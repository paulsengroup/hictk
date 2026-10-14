// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <type_traits>

namespace hictk::cooler {

template <typename PixelT, typename>
inline Attributes Attributes::init(std::uint32_t bin_size_) {
  Attributes attrs{};
  attrs.bin_size = bin_size_;
  if constexpr (std::is_floating_point_v<PixelT>) {
    attrs.sum = 0.0;
    attrs.cis = 0.0;
  }
  if constexpr (std::is_integral_v<PixelT>) {
    attrs.sum = std::int64_t{0};
    attrs.cis = std::int64_t{0};
  }
  return attrs;
}

}  // namespace hictk::cooler
