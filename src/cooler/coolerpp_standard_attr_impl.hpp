// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <type_traits>

namespace coolerpp {

template <typename PixelT, typename>
inline StandardAttributes StandardAttributes::init(std::uint32_t bin_size_) {
  StandardAttributes attrs{};
  attrs.bin_size = bin_size_;
  if constexpr (std::is_floating_point_v<PixelT>) {
    attrs.sum = 0.0;
    attrs.cis = 0.0;
  }
  if constexpr (std::is_integral_v<PixelT>) {
    attrs.sum = std::int64_t(0);
    attrs.cis = std::int64_t(0);
  }
  return attrs;
}

inline StandardAttributes StandardAttributes::init_empty() noexcept {
  StandardAttributes attrs{};

  attrs.bin_type.reset();
  attrs.creation_date.reset();
  attrs.format_url.reset();
  attrs.generated_by.reset();
  attrs.assembly.reset();
  attrs.nbins.reset();
  attrs.nchroms.reset();
  attrs.metadata.reset();
  attrs.storage_mode.reset();
  attrs.sum.reset();
  attrs.cis.reset();

  return attrs;
}

inline bool StandardAttributes::operator==(
    const coolerpp::StandardAttributes& other) const noexcept {
  // clang-format off
  return bin_size == other.bin_size &&
         bin_type == other.bin_type &&
         format == other.format &&
         format_version == other.format_version &&
         storage_mode == other.storage_mode &&
         creation_date == other.creation_date &&
         generated_by == other.generated_by &&
         assembly == other.assembly &&
         metadata == other.metadata &&
         format_url == other.format_url &&
         nbins == other.nbins &&
         nchroms == other.nchroms &&
         nnz == other.nnz &&
         sum == other.sum &&
         cis == other.cis;
  // clang-format on
}

inline bool StandardAttributes::operator!=(
    const coolerpp::StandardAttributes& other) const noexcept {
  return !(*this == other);
}
}  // namespace coolerpp
