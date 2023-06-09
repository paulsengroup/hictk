// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <functional>
#include <iterator>

namespace hictk::internal {

constexpr HiCHeader::operator bool() const noexcept { return masterIndexOffset >= 0; }

inline bool HiCHeader::operator==(const HiCHeader &other) const noexcept {
  return url == other.url && masterIndexOffset == other.masterIndexOffset;
}

inline bool HiCHeader::operator!=(const HiCHeader &other) const noexcept {
  return !(*this == other);
}

}  // namespace hictk::internal

template <>
struct std::hash<hictk::internal::HiCHeader> {
  inline std::size_t operator()(hictk::internal::HiCHeader const &h) const noexcept {
    return hictk::internal::hash_combine(0, h.url, h.masterIndexOffset);
  }
};
