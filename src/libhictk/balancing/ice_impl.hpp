// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
#include <type_traits>

#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::balancing {

template <typename PixelIt>
inline ICE::ICE(PixelIt first_pixel, PixelIt last_pixel, std::size_t num_rows, double tol,
                std::size_t max_iters, std::size_t num_masked_diags, std::size_t min_nnz,
                [[maybe_unused]] double min_count)
    : _biases(num_rows, 1.0) {
  auto [bin1_ids, bin2_ids, counts] =
      construct_sparse_matrix(first_pixel, last_pixel, num_masked_diags);
  std::vector<double> margs(_biases.size());

  if (min_nnz != 0) {
    filter_rows_by_nnz(bin1_ids, bin2_ids, counts, _biases, min_nnz, margs);
  }

  // TODO mad-max filter

  for (std::size_t i = 0; i < max_iters; ++i) {
    const auto res = inner_loop(bin1_ids, bin2_ids, counts, _biases, margs);
    _variance = res.variance;
    _scale = res.scale;
    if (res.variance < tol) {
      break;
    }
  }
}

template <typename PixelIt>
inline std::tuple<std::vector<std::size_t>, std::vector<std::size_t>, std::vector<double>>
ICE::construct_sparse_matrix(PixelIt first_pixel, PixelIt last_pixel,
                             std::size_t num_masked_diags) {
  std::vector<std::size_t> bin1_ids{};
  std::vector<std::size_t> bin2_ids{};
  std::vector<double> counts{};
  std::for_each(first_pixel, last_pixel, [&](const auto& p) {
    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
      bin1_ids.push_back(p.bin1_id);
      bin2_ids.push_back(p.bin2_id);
      counts.push_back(p.count);
    }
  });

  bin1_ids.shrink_to_fit();
  bin2_ids.shrink_to_fit();
  counts.shrink_to_fit();
  return std::make_tuple(bin1_ids, bin2_ids, counts);
}

inline void ICE::times_outer_product(const std::vector<std::size_t>& bin1_ids,
                                     const std::vector<std::size_t>& bin2_ids,
                                     std::vector<double>& counts,
                                     const std::vector<double>& biases) {
  assert(bin1_ids.size() == counts.size());
  assert(bin2_ids.size() == counts.size());
  for (std::size_t i = 0; i < counts.size(); ++i) {
    const auto i1 = bin1_ids[i];
    const auto i2 = bin2_ids[i];
    counts[i] *= biases[i1] * biases[i2];
  }
}

inline void ICE::marginalize(const std::vector<std::size_t>& bin1_ids,
                             const std::vector<std::size_t>& bin2_ids, std::vector<double>& counts,
                             std::vector<double>& marg) {
  std::fill(marg.begin(), marg.end(), 0);

  for (std::size_t i = 0; i < counts.size(); ++i) {
    const auto i1 = bin1_ids[i];
    const auto i2 = bin2_ids[i];

    marg[i1] += counts[i];
    marg[i2] += counts[i];
  }
}

inline void ICE::filter_rows_by_nnz(const std::vector<std::size_t>& bin1_ids,
                                    const std::vector<std::size_t>& bin2_ids,
                                    std::vector<double> counts, std::vector<double>& biases,
                                    std::size_t min_nnz, std::vector<double>& marg_buff) {
  std::transform(counts.begin(), counts.end(), counts.begin(), [](const auto n) { return n != 0; });
  marginalize(bin1_ids, bin2_ids, counts, marg_buff);
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (marg_buff[i] < static_cast<double>(min_nnz)) {
      biases[i] = 0;
    }
  }
}

inline auto ICE::inner_loop(const std::vector<std::size_t>& bin1_ids,
                            const std::vector<std::size_t>& bin2_ids, std::vector<double> counts,
                            std::vector<double>& biases, std::vector<double>& marg_buffer)
    -> Result {
  times_outer_product(bin1_ids, bin2_ids, counts, biases);

  marginalize(bin1_ids, bin2_ids, counts, marg_buffer);

  double marg_sum = 0.0;
  std::size_t nnz_marg{};
  for (const auto& n : marg_buffer) {
    marg_sum += n;
    nnz_marg += n != 0;
  }

  if (nnz_marg == 0) {
    std::fill(biases.begin(), biases.end(), std::numeric_limits<double>::quiet_NaN());
    return {std::numeric_limits<double>::quiet_NaN(), 0.0};
  }

  const auto avg_nzmarg = (marg_sum / static_cast<double>(nnz_marg));
  for (std::size_t i = 0; i < biases.size(); ++i) {
    const auto n = marg_buffer[i] / avg_nzmarg;
    if (n != 0) {
      biases[i] /= n;
    }
  }

  double ssq_nzmarg = 0.0;
  for (const auto n : marg_buffer) {
    if (n != 0) {
      ssq_nzmarg += std::pow(n - avg_nzmarg, 2);
    }
  }
  const auto var_nzmarg = ssq_nzmarg / static_cast<double>(nnz_marg - 1);

  return {avg_nzmarg, var_nzmarg};
}

inline std::vector<double> ICE::get_weights(bool rescale) const {
  std::vector<double> biases(_biases.size());
  const auto scale = rescale ? std::sqrt(_scale) : 1.0;
  std::transform(_biases.begin(), _biases.end(), biases.begin(), [&](const auto n) {
    return n == 0 ? std::numeric_limits<double>::quiet_NaN() : n / scale;
  });
  return biases;
}

inline double ICE::scale() const noexcept { return _scale; }
inline double ICE::variance() const noexcept { return _variance; }

}  // namespace hictk::balancing
