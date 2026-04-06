// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/expected_values_aggregator.hpp"

#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"

namespace hictk {

ExpectedValuesAggregator::ExpectedValuesAggregator(std::shared_ptr<const BinTable> bins)
    : _bins(std::move(bins)) {
  SPDLOG_INFO(FMT_STRING("[{} bp] initializing expected value vector"), _bins->resolution());
  std::uint32_t max_length = 0;
  for (const auto &chrom : chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }

    max_length = std::max(max_length, chrom.size());

    _num_bins_gw += chrom.size();
  }

  const auto bin_size = _bins->resolution();
  // round down to mimic HiCTools' behavior
  const auto max_n_bins = max_length / bin_size;
  _possible_distances.resize(max_n_bins, 0.0);
  _actual_distances.resize(max_n_bins, 0.0);
  _weights.resize(max_n_bins, 0.0);
}

void ExpectedValuesAggregator::compute_density() {
  SPDLOG_INFO(FMT_STRING("[{} bp] computing expected vector density"), _bins->resolution());
  init_possible_distances();
  compute_density_cis();
  compute_density_trans();
}

const std::vector<double> &ExpectedValuesAggregator::weights() const noexcept { return _weights; }

std::vector<double> ExpectedValuesAggregator::weights(const Chromosome &chrom, bool rescale) const {
  // We round down to mach HiCTools behavior
  const auto num_bins = chrom.size() / _bins->resolution();

  std::vector<double> w{_weights.begin(), _weights.begin() + static_cast<std::ptrdiff_t>(num_bins)};
  if (!rescale) {
    return w;
  }

  const auto sf = scaling_factor(chrom);
  std::transform(w.begin(), w.end(), w.begin(), [&](const auto n) { return n / sf; });

  return w;
}

double ExpectedValuesAggregator::scaling_factor(const Chromosome &chrom) const {
  return _scaling_factors.at(chrom);
}

const phmap::btree_map<Chromosome, double> &ExpectedValuesAggregator::scaling_factors()
    const noexcept {
  return _scaling_factors;
}

void ExpectedValuesAggregator::init_possible_distances() {
  const auto bin_size = _bins->resolution();

  for (const auto &[chrom, _] : _cis_sum) {
    if (chrom.is_all()) {
      continue;
    }
    const auto n_bins = chrom.size() / bin_size;
    for (std::uint32_t i = 0; i < n_bins; ++i) {
      _possible_distances[i] += n_bins - i;
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
void ExpectedValuesAggregator::compute_density_cis() {
  // Re-implementation of the algorithm used by HiCTools:
  // https://github.com/aidenlab/HiCTools/blob/6b2fab8e78685deae199c33bbb167dcab1dbfbb3/src/hic/tools/utils/original/ExpectedValueCalculation.java#L184

  if (_actual_distances.empty()) {
    return;
  }

  auto num_sum = _actual_distances.front();
  auto den_sum = _possible_distances.front();
  std::size_t bound1 = 0;
  std::size_t bound2 = 0;

  const auto shot_noise_minimum = 400.0;
  const auto max_num_bins = _actual_distances.size();

  assert(_weights.size() == max_num_bins);
  for (std::size_t ii = 0; ii < max_num_bins; ii++) {
    if (num_sum < shot_noise_minimum) {
      while (num_sum < shot_noise_minimum && ++bound2 < max_num_bins) {
        num_sum += _actual_distances[bound2];
        den_sum += _possible_distances[bound2];
      }
    } else if (num_sum >= shot_noise_minimum && bound2 > bound1) {
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

  for (const auto &[chrom, _] : _cis_sum) {
    if (chrom.is_all()) {
      continue;
    }
    auto num_chrom_bins = chrom.size() / _bins->resolution();
    auto expected_count = 0.0;
    for (std::size_t n = 0; n < num_chrom_bins; n++) {
      if (n < max_num_bins) {
        const double v = _weights[n];
        expected_count += (static_cast<double>(num_chrom_bins) - static_cast<double>(n)) * v;
      }
    }

    const double observed_count = _cis_sum.at(chrom);
    const double f = expected_count / observed_count;
    _scaling_factors.emplace(chrom, f);
  }
}

void ExpectedValuesAggregator::compute_density_trans() {
  for (auto &[k, v] : _trans_sum) {
    // We round-down to match HiCTools behavior
    const auto num_bins1 = k.first.size() / _bins->resolution();
    const auto num_bins2 = k.second.size() / _bins->resolution();
    const auto num_pixels = num_bins1 * num_bins2;
    v = num_pixels != 0 ? v / static_cast<double>(num_pixels) : 0.0;
  }
}

double ExpectedValuesAggregator::at(const Chromosome &chrom) const { return _cis_sum.at(chrom); }

double ExpectedValuesAggregator::at(const Chromosome &chrom1, const Chromosome &chrom2) const {
  return _trans_sum.at(std::make_pair(chrom1, chrom2));
}

double &ExpectedValuesAggregator::at(const Chromosome &chrom) {
  auto [it, _] = _cis_sum.try_emplace(chrom, 0.0);
  return it->second;
}

double &ExpectedValuesAggregator::at(const Chromosome &chrom1, const Chromosome &chrom2) {
  auto [it, _] = _trans_sum.try_emplace(std::make_pair(chrom1, chrom2), 0.0);
  return it->second;
}

const Reference &ExpectedValuesAggregator::chromosomes() const noexcept {
  assert(_bins);
  return _bins->chromosomes();
}

}  // namespace hictk
