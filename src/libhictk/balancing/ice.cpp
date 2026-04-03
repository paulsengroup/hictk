// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/balancing/ice.hpp"

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <mutex>
#include <nonstd/span.hpp>
#include <utility>
#include <vector>

#include "hictk/balancing/weights.hpp"

namespace hictk::balancing {

void ICE::min_count_filtering(nonstd::span<double> biases, std::size_t min_count,
                              nonstd::span<const double> marg) {
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (marg[i] < static_cast<double>(min_count)) {
      biases[i] = 0;
    }
  }
}

void ICE::mad_max_filtering(nonstd::span<const std::uint64_t> chrom_offsets,
                            nonstd::span<double> biases, nonstd::span<double> marg,
                            double mad_max) {
  auto median = [](auto v) {
    assert(!v.empty());

    const auto size = static_cast<std::ptrdiff_t>(v.size());
    auto first = v.begin();
    auto mid = first + (size / 2);
    auto last = v.end();

    std::nth_element(first, mid, last);

    if (size % 2 != 0) {
      return *mid;
    }

    const auto n1 = *mid;
    std::nth_element(first, --mid, last);
    const auto n2 = *mid;

    return (n1 + n2) / 2;
  };

  auto mad = [&](const auto& vin) {
    const auto median_ = median(vin);
    auto vout = vin;

    std::transform(vout.begin(), vout.end(), vout.begin(),
                   [&](const auto n) { return std::abs(n - median_); });

    return median(vout);
  };

  assert(chrom_offsets.size() > 1);
  std::vector<double> cmarg{};
  for (std::size_t i = 1; i < chrom_offsets.size(); ++i) {
    const auto i0 = static_cast<std::ptrdiff_t>(chrom_offsets[i - 1] - chrom_offsets.front());
    const auto i1 = static_cast<std::ptrdiff_t>(chrom_offsets[i] - chrom_offsets.front());

    cmarg.clear();  // NOLINTNEXTLINE(*-pointer-arithmetic)
    std::copy_if(marg.begin() + i0, marg.begin() + i1, std::back_inserter(cmarg),
                 [](const auto n) { return n > 0; });

    if (!cmarg.empty()) {
      const auto median_ = median(cmarg);  // NOLINTNEXTLINE(*-pointer-arithmetic)
      std::transform(marg.begin() + i0, marg.begin() + i1, marg.begin() + i0,
                     [&](const auto n) { return n / median_; });
    }
  }

  std::vector<double> log_nz_marg{};
  for (const auto n : marg) {
    if (n > 0) {
      log_nz_marg.push_back(std::log(n));
    }
  }

  if (log_nz_marg.empty()) {
    return;
  }

  const auto median_log_nz_marg = median(log_nz_marg);
  const auto dev_log_nz_marg = mad(log_nz_marg);

  const auto cutoff = std::exp(median_log_nz_marg - (mad_max * dev_log_nz_marg));

  for (std::size_t i = 0; i < marg.size(); ++i) {
    if (marg[i] < cutoff) {
      biases[i] = 0.0;
    }
  }
}

std::pair<double, std::size_t> ICE::aggregate_marg(nonstd::span<const double> marg,
                                                   BS::light_thread_pool* tpool) {
  double marg_sum = 0.0;
  std::size_t nnz_marg{};

  std::mutex mtx{};

  auto aggregate_marg_impl = [&](std::size_t istart, std::size_t iend) {
    double marg_sum_ = 0.0;
    std::size_t nnz_marg_{};

    for (auto i = istart; i < iend; ++i) {
      marg_sum_ += marg[i];
      nnz_marg_ += static_cast<std::size_t>(marg[i] != 0);
    }

    [[maybe_unused]] const std::scoped_lock lck(mtx);
    marg_sum += marg_sum_;
    nnz_marg += nnz_marg_;
  };

  if (!process_in_parallel(marg, tpool)) {
    aggregate_marg_impl(0, marg.size());
    return std::make_pair(marg_sum, nnz_marg);
  }

  tpool->detach_blocks(std::size_t{0}, marg.size(), aggregate_marg_impl);
  tpool->wait();

  return std::make_pair(marg_sum, nnz_marg);
}

void ICE::update_biases(nonstd::span<const double> marg, nonstd::span<double> biases,
                        double avg_nzmarg, BS::light_thread_pool* tpool) {
  auto update_biases_impl = [&](std::size_t istart, std::size_t iend) {
    for (auto i = istart; i < iend; ++i) {
      const auto n = marg[i] / avg_nzmarg;
      if (n != 0) {
        biases[i] /= n;
      }
    }
  };

  if (!process_in_parallel(marg, tpool)) {
    update_biases_impl(0, marg.size());
    return;
  }

  tpool->detach_blocks(std::size_t{0}, marg.size(), update_biases_impl);
  tpool->wait();
}

double ICE::compute_ssq_nzmarg(nonstd::span<const double> marg, double avg_nzmarg,
                               BS::light_thread_pool* tpool) {
  std::mutex mtx{};
  double ssq_nzmarg = 0;

  auto compute_ssq_nzmarg_impl = [&](std::size_t istart, std::size_t iend) {
    double ssq_nzmarg_ = 0.0;
    for (auto i = istart; i < iend; ++i) {
      const auto& n = marg[i];
      if (n != 0) {
        ssq_nzmarg_ += std::pow(n - avg_nzmarg, 2);
      }
    }
    [[maybe_unused]] const std::scoped_lock lck(mtx);
    ssq_nzmarg += ssq_nzmarg_;
  };

  if (!process_in_parallel(marg, tpool)) {
    compute_ssq_nzmarg_impl(0, marg.size());
    return ssq_nzmarg;
  }

  tpool->detach_blocks(std::size_t{0}, marg.size(), compute_ssq_nzmarg_impl);
  tpool->wait();
  return ssq_nzmarg;
}

std::vector<double> ICE::compute_weights_from_chromosome_sizes(
    const BinTable& bins, nonstd::span<std::uint64_t> chrom_bin_offsets) {
  std::vector<double> weights(bins.size());
  for (std::uint32_t i = 1; i < chrom_bin_offsets.size(); ++i) {
    const auto& chrom = bins.chromosomes().at(i - 1);
    if (chrom.is_all()) {
      continue;
    }
    const auto i0 = chrom_bin_offsets[i - 1];
    const auto i1 = chrom_bin_offsets[i];

    const auto nbins = static_cast<double>(bins.size());
    const auto cnbins =
        std::ceil(static_cast<double>(chrom.size()) / static_cast<double>(bins.resolution()));

    for (std::size_t j = i0; j < i1; ++j) {
      weights[j] = 1.0 / (1.0 - cnbins / nbins);
    }
  }
  return weights;
}

balancing::Weights ICE::get_weights(bool rescale) const {
  if (!rescale) {
    return {_biases, balancing::Weights::Type::MULTIPLICATIVE};
  }

  std::vector<double> biases(_biases.size());
  if (_scale.size() == 1) {
    const auto scale = std::sqrt(_scale[0]);
    std::transform(_biases.begin(), _biases.end(), biases.begin(), [&](const auto n) {
      return n == 0 ? std::numeric_limits<double>::quiet_NaN() : n / scale;
    });
  } else {
    for (std::size_t i = 1; i < _chrom_offsets.size(); ++i) {
      const auto i0 = static_cast<std::ptrdiff_t>(_chrom_offsets[i - 1]);
      const auto i1 = static_cast<std::ptrdiff_t>(_chrom_offsets[i]);
      const auto scale = std::sqrt(_scale[i - 1]);
      std::transform(_biases.begin() + i0, _biases.begin() + i1, biases.begin() + i0,
                     [&](const auto n) {
                       return n == 0 ? std::numeric_limits<double>::quiet_NaN() : n / scale;
                     });
    }
  }
  return {biases, balancing::Weights::Type::MULTIPLICATIVE};
}

std::vector<double> ICE::scale() const noexcept { return _scale; }
std::vector<double> ICE::variance() const noexcept { return _variance; }

}  // namespace hictk::balancing
