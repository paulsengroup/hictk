// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

inline ExpectedValuesAggregator::ExpectedValuesAggregator(std::shared_ptr<const BinTable> bins)
    : _bins(std::move(bins)) {
  SPDLOG_INFO(FMT_STRING("[{} bp] initializing expected value vector"), _bins->bin_size());
  std::uint32_t max_length = 0;
  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }

    max_length = std::max(max_length, chrom1.size());
    _cis_sum.emplace(chrom1, 0.0);

    for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);
      _trans_sum.emplace(std::make_pair(chrom1, chrom2), 0.0);
    }

    for (const auto &chrom : chromosomes()) {
      if (chrom.is_all()) {
        continue;
      }
      _num_bins_gw += chrom.size();
    }
  }

  const auto bin_size = _bins->bin_size();
  const auto max_n_bins = (max_length + bin_size - 1) / bin_size;
  _possible_distances.resize(max_n_bins, 0.0);
  _actual_distances.resize(max_n_bins, 0.0);

  for (const auto &chrom : chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto n_bins = chrom.size() / bin_size;
    for (std::uint32_t i = 0; i < n_bins; ++i) {
      _possible_distances[i] += n_bins - i;
    }
  }
}

inline void ExpectedValuesAggregator::add(const ThinPixel<float> &p) {
  add(Pixel<float>{*_bins, p});
}

inline void ExpectedValuesAggregator::add(const Pixel<float> &p) {
  if (std::isnan(p.count)) {
    return;
  }

  const auto &chrom1 = p.coords.bin1.chrom();
  const auto &chrom2 = p.coords.bin2.chrom();

  if (p.coords.is_intra()) {
    at(chrom1) += static_cast<double>(p.count);
    const auto i = p.coords.bin2.id() - p.coords.bin1.id();
    _actual_distances[i] += static_cast<double>(p.count);
  } else {
    at(chrom1, chrom2) += static_cast<double>(p.count);
  }
}

inline void ExpectedValuesAggregator::compute_density() {
  SPDLOG_INFO(FMT_STRING("[{} bp] computing expected vector density"), _bins->bin_size());
  compute_density_cis();
  compute_density_trans();
}

inline const std::vector<double> &ExpectedValuesAggregator::weights() const noexcept {
  return _weights;
}

inline double ExpectedValuesAggregator::scaling_factor(const Chromosome &chrom) const {
  return _scaling_factors.at(chrom);
}

inline const phmap::btree_map<Chromosome, double> &ExpectedValuesAggregator::scaling_factors()
    const noexcept {
  return _scaling_factors;
}

inline void ExpectedValuesAggregator::compute_density_cis() {
  // Re-implementation of the algorithm used by HiCTools:
  // https://github.com/aidenlab/HiCTools/blob/6b2fab8e78685deae199c33bbb167dcab1dbfbb3/src/hic/tools/utils/original/ExpectedValueCalculation.java#L184

  auto num_sum = _actual_distances.front();
  auto den_sum = _possible_distances.front();
  std::size_t bound1 = 0;
  std::size_t bound2 = 0;

  const auto shot_noise_minimum = 400.0;
  const auto max_num_bins = _actual_distances.size();

  _weights.resize(max_num_bins);
  std::fill(_weights.begin(), _weights.end(), 0.0);

  for (std::size_t ii = 0; ii < max_num_bins; ii++) {
    if (num_sum < shot_noise_minimum) {
      while (num_sum < shot_noise_minimum && ++bound2 < max_num_bins) {
        num_sum += _actual_distances[bound2];
        den_sum += _possible_distances[bound2];
      }
    } else if (num_sum >= shot_noise_minimum && bound2 - bound1 > 0) {
      while (bound2 > bound1 && bound2 < _num_bins_gw && bound1 < _num_bins_gw &&
             num_sum - _actual_distances[bound1] - _actual_distances[bound2] >=
                 shot_noise_minimum) {
        num_sum = num_sum - _actual_distances[bound1] - _actual_distances[bound2];
        den_sum = den_sum - _possible_distances[bound1] - _possible_distances[bound2];
        bound1++;
        bound2--;
      }
    }
    _weights[ii] = num_sum / den_sum;
    if (bound2 + 2 < max_num_bins) {
      num_sum += _actual_distances[bound2 + 1] + _actual_distances[bound2 + 2];
      den_sum += _possible_distances[bound2 + 1] + _possible_distances[bound2 + 2];
      bound2 += 2;
    } else if (bound2 + 1 < max_num_bins) {
      num_sum += _actual_distances[bound2 + 1];
      den_sum += _possible_distances[bound2 + 1];
      bound2++;
    }
  }

  for (const auto &chrom : chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    auto num_chrom_bins = chrom.size() / _bins->bin_size();
    auto expected_count = 0.0;
    for (std::size_t n = 0; n < num_chrom_bins; n++) {
      if (n < max_num_bins) {
        double v = _weights[n];
        expected_count += (double(num_chrom_bins) - double(n)) * v;
      }
    }

    double observed_count = _cis_sum.at(chrom);
    double f = expected_count / observed_count;
    _scaling_factors.emplace(chrom, f);
  }
}

inline void ExpectedValuesAggregator::compute_density_trans() {
  for (auto &[k, v] : _trans_sum) {
    // We round-down to match HiCTools behavior
    const auto num_bins1 = k.first.size() / _bins->bin_size();
    const auto num_bins2 = k.second.size() / _bins->bin_size();
    const auto num_pixels = num_bins1 * num_bins2;
    v = num_pixels != 0 ? v / static_cast<double>(num_pixels) : 0.0;
  }
}

inline double ExpectedValuesAggregator::at(const Chromosome &chrom) const {
  return _cis_sum.at(chrom);
}

inline double ExpectedValuesAggregator::at(const Chromosome &chrom1,
                                           const Chromosome &chrom2) const {
  return _trans_sum.at(std::make_pair(chrom1, chrom2));
}

inline double &ExpectedValuesAggregator::at(const Chromosome &chrom) { return _cis_sum.at(chrom); }

inline double &ExpectedValuesAggregator::at(const Chromosome &chrom1, const Chromosome &chrom2) {
  return _trans_sum.at(std::make_pair(chrom1, chrom2));
}

inline const Reference &ExpectedValuesAggregator::chromosomes() const noexcept {
  assert(_bins);
  return _bins->chromosomes();
}

}  // namespace hictk::hic::internal
