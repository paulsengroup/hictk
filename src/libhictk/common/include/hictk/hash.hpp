// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <functional>

namespace hictk::internal {
// Adapted from:
// https://www.boost.org/doc/libs/1_37_0/doc/html/hash/reference.html#boost.hash_combine

template <typename T>
[[nodiscard]] inline std::size_t hash_combine(std::size_t seed, const T &v) {
  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
  return seed;
}
template <typename T, typename... Args>
[[nodiscard]] inline std::size_t hash_combine(std::size_t seed, const T &v, const Args &...args) {
  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
  return hash_combine(seed, args...);
}
}  // namespace hictk::internal
