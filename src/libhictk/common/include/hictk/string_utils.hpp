// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <string_view>

namespace hictk::internal {

/// Checks whether a string starts with the given prefix
[[nodiscard]] constexpr bool starts_with(std::string_view s, std::string_view prefix) {
  if (s.size() < prefix.size()) {
    return false;
  }
  for (std::size_t i = 0; i < prefix.size(); ++i) {
    if (s[i] != prefix[i]) {
      return false;
    }
  }
  return true;
}

/// Checks whether a string ends with the given suffix
[[nodiscard]] constexpr bool ends_with(std::string_view s, std::string_view suffix) {
  if (s.size() < suffix.size()) {
    return false;
  }
  for (std::size_t i = 1; i <= suffix.size(); ++i) {
    const auto i1 = s.size() - i;
    const auto i2 = suffix.size() - i;
    if (s[i1] != suffix[i2]) {
      return false;
    }
  }
  return true;
}

}  // namespace hictk::internal
