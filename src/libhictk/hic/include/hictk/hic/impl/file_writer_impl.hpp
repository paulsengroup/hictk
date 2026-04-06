// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iosfwd>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "hictk/type_traits.hpp"

namespace hictk::hic::internal {

template <typename I1, typename I2>
HiCSectionOffsets::HiCSectionOffsets(I1 start_, I2 size_)
    : _position(conditional_static_cast<std::streampos>(start_)),
      _size(conditional_static_cast<std::size_t>(size_)) {
  static_assert(std::is_integral_v<I1> || is_specialization_v<I1, std::fpos>);
  static_assert(std::is_integral_v<I2> || is_specialization_v<I2, std::fpos>);

  if constexpr (std::is_signed_v<I1> || is_specialization_v<I1, std::fpos>) {
    if (start_ < 0) {
      throw std::logic_error(fmt::format(
          FMT_STRING("start position for HiCSectionOffset cannot be negative, found {}"),
          static_cast<std::int64_t>(start_)));
    }
  }
  if constexpr (std::is_signed_v<I2> || is_specialization_v<I2, std::fpos>) {
    assert(size_ >= 0);  // NOLINT
    if (size_ < 0) {
      throw std::logic_error(
          fmt::format(FMT_STRING("size given to HiCSectionOffset cannot be negative, found {}"),
                      static_cast<std::int64_t>(start_)));
    }
  }
}

template <typename PixelIt, typename>
void HiCFileWriter::add_pixels(std::uint32_t resolution, PixelIt first_pixel, PixelIt last_pixel,
                               bool validate) {
  try {
    assert(_tpool);
    _block_mappers.at(resolution)
        .append_pixels(std::move(first_pixel), std::move(last_pixel), validate, *_tpool);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while adding pixels for resolution {} to file \"{}\": {}"),
        resolution, path(), e.what()));
  }
}

}  // namespace hictk::hic::internal
