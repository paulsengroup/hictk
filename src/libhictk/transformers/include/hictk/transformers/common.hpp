// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <type_traits>

#include "hictk/pixel.hpp"
#include "hictk/suppress_warnings.hpp"

namespace hictk::transformers {
enum class QuerySpan : std::uint_fast8_t { lower_triangle, upper_triangle, full };

namespace internal {

template <typename T, typename = std::void_t<>>
inline constexpr bool has_coord1_member_fx = false;

template <typename T>
inline constexpr bool has_coord1_member_fx<T, std::void_t<decltype(std::declval<T>().coord1())>> =
    true;

template <typename N, typename PixelSelector, typename MatrixT, typename SetterOp>
inline void fill_matrix(const PixelSelector& sel, MatrixT& buffer, std::int64_t num_rows,
                        std::int64_t num_cols, std::int64_t offset1, std::int64_t offset2,
                        bool populate_lower_triangle, bool populate_upper_triangle,
                        SetterOp matrix_setter) {
  assert(populate_lower_triangle || populate_upper_triangle);

  std::for_each(sel.template begin<N>(), sel.template end<N>(), [&](const ThinPixel<N>& p) {
    const auto i1 = static_cast<std::int64_t>(p.bin1_id) - offset1;
    const auto i2 = static_cast<std::int64_t>(p.bin2_id) - offset2;
    bool inserted = false;
    if (populate_upper_triangle) {
      if (i1 >= 0 && i1 < num_rows && i2 >= 0 && i2 < num_cols) {
        HICTK_DISABLE_WARNING_PUSH
        HICTK_DISABLE_WARNING_CONVERSION
        matrix_setter(buffer, i1, i2, p.count);
        HICTK_DISABLE_WARNING_POP
        inserted = true;
      }
    }

    if (populate_lower_triangle) {
      const auto i3 = static_cast<std::int64_t>(p.bin2_id) - offset1;
      const auto i4 = static_cast<std::int64_t>(p.bin1_id) - offset2;

      if (inserted && i1 == i3 && i2 == i4) {
        return;
      }

      if (i3 >= 0 && i3 < num_rows && i4 >= 0 && i4 < num_cols) {
        HICTK_DISABLE_WARNING_PUSH
        HICTK_DISABLE_WARNING_CONVERSION
        matrix_setter(buffer, i3, i4, p.count);
        HICTK_DISABLE_WARNING_POP
      }
    }
  });
}

}  // namespace internal
}  // namespace hictk::transformers
