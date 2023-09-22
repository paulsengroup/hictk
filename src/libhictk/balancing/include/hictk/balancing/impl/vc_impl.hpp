// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>
#include <iostream>
#include <iterator>
#include <type_traits>

#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::balancing {

template <typename N>
template <typename PixelIt>
inline VC<N>::VC(PixelIt first_pixel, PixelIt last_pixel, std::size_t num_rows,
                 std::size_t bin_id_offset) {
  if constexpr (std::is_floating_point_v<N>) {
    _rowsum = std::vector<double>(num_rows, 0);
    _sum = 0.0;
  } else {
    _rowsum = std::vector<std::int64_t>(num_rows, 0);
    _sum = std::int64_t(0);
  }

  // Compute rowsum and matrix sum
  std::visit(
      [&](auto& sum) {
        using T = remove_cvref_t<decltype(sum)>;
        auto& rowsum = std::get<std::vector<T>>(_rowsum);
        std::for_each(first_pixel, last_pixel, [&](const ThinPixel<N>& p) {
          if constexpr (std::is_floating_point_v<N>) {
            if (std::isnan(p.count)) {
              return;
            }
          }
          const auto bin1_id = p.bin1_id - bin_id_offset;
          const auto bin2_id = p.bin2_id - bin_id_offset;
          const auto count = conditional_static_cast<T>(p.count);

          rowsum[bin1_id] += count;

          if (bin1_id != bin2_id) {
            rowsum[bin2_id] += count;
          }
        });
      },
      _sum);

  std::visit(
      [&](auto& sum) {
        using T = remove_cvref_t<decltype(sum)>;
        const auto& rowsum = std::get<std::vector<T>>(_rowsum);
        std::for_each(first_pixel, last_pixel, [&](const ThinPixel<N>& p) {
          if constexpr (std::is_floating_point_v<N>) {
            if (std::isnan(p.count)) {
              return;
            }
          }
          const auto bin1_id = p.bin1_id - bin_id_offset;
          const auto bin2_id = p.bin2_id - bin_id_offset;

          const auto rs1 = conditional_static_cast<double>(rowsum[bin1_id]);
          const auto rs2 = conditional_static_cast<double>(rowsum[bin2_id]);
          if (rs1 == 0 || rs2 == 0) {
            return;
          }

          const auto count = conditional_static_cast<T>(bin1_id == bin2_id ? p.count : 2 * p.count);
          sum += count;
          _norm_sum += conditional_static_cast<double>(count) / (rs1 * rs2);
        });
      },
      _sum);
}

template <typename N>
inline std::vector<double> VC<N>::get_weights() const {
  std::vector<double> weights;

  const auto scaling_factor = std::visit(
      [&](const auto& sum) { return std::sqrt(_norm_sum / conditional_static_cast<double>(sum)); },
      _sum);

  std::visit(
      [&](const auto& rowsum) {
        weights.reserve(rowsum.size());

        for (const auto rs : rowsum) {
          weights.push_back(conditional_static_cast<double>(rs) * scaling_factor);
        }
      },
      _rowsum);

  return weights;
}

}  // namespace hictk::balancing
