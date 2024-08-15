// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstddef>
#include <string>
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

/// Checks whether a string ends with the given prefix
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

/// Remove the given prefix from a string
[[nodiscard]] constexpr std::string_view remove_prefix(std::string_view s,
                                                       std::string_view suffix) {
  if (starts_with(s, suffix)) {
    s.remove_prefix(suffix.size());
  }
  return s;  // NOLINT
}

/// Remove the given suffix from a string
[[nodiscard]] constexpr std::string_view remove_suffix(std::string_view s,
                                                       std::string_view suffix) {
  if (ends_with(s, suffix)) {
    s.remove_suffix(suffix.size());
  }
  return s;  // NOLINT
}

/// Strip the first (i.e. outer) quote pair
[[nodiscard]] constexpr std::string_view strip_first_quote_pair(std::string_view s,
                                                                char quote_symbol = '"') {
  if (s.front() == quote_symbol && s.back() == quote_symbol && s.size() > 1) {
    return s.substr(1, s.size() - 2);
  }

  return s;  // NOLINT
}

/// Escape string
[[nodiscard]] inline std::string escape_str(std::string_view s) {
  const auto str = fmt::format(FMT_STRING("{:?}"), s);
  return str.substr(1, str.size() - 2);
}

[[nodiscard]] inline std::string str_replace(std::string_view s, std::string_view old_str,
                                             std::string_view new_str) {
  std::string ss{s};

  std::size_t pos = 0;
  while ((pos = ss.find(old_str, pos)) != std::string::npos) {
    ss.replace(pos, old_str.length(), new_str);
    pos += new_str.length();
  }
  return ss;
}

}  // namespace hictk::internal
