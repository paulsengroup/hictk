// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fast_float/fast_float.h>  // for from_chars (fp)
#include <fmt/format.h>             // for compile_string_to_view, FMT_STRING

#include <charconv>      // for from_chars (int)
#include <limits>        // for numeric_limits
#include <stdexcept>     // for runtime_error, logic_error
#include <string>        // for string
#include <string_view>   // for string_view
#include <system_error>  // for errc, make_error_code, errc::invalid_argument, errc:...
#include <type_traits>   // for is_arithmetic, is_integral, is_unsigned
#include <vector>        // for vector

#include "hictk/type_pretty_printer.hpp"

namespace hictk::internal {

template <typename N>
inline void throw_except_from_errc(std::string_view tok, std::size_t idx,
                                   [[maybe_unused]] const N &field, const char *c, std::errc e) {
  static_assert(std::is_arithmetic<N>());
  std::string base_error;
  if (idx != (std::numeric_limits<std::size_t>::max)()) {
    base_error = fmt::format(FMT_STRING("Unable to convert field {} (\"{}\") to a "), idx, tok);
  } else {
    base_error = fmt::format(FMT_STRING("Unable to convert field \"{}\" to"), tok);
  }
  if (std::is_integral<N>()) {
    if (std::is_unsigned<N>()) {
      base_error += " a positive integral number";
    } else {
      base_error += " an integral number";
    }
  } else {
    base_error += " a real number";
  }
  if (e == std::errc::invalid_argument) {
    if (c != nullptr) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("{}. Reason: found an invalid character \"{}\""), base_error, *c));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("{}. Reason: found an invalid character"), base_error));
  }
  if (e == std::errc::result_out_of_range) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("{}. Reason: number {} is outside the range of representable numbers [{}, {}]."),
        base_error, tok, (std::numeric_limits<N>::min)(), (std::numeric_limits<N>::max)()));
  }

  throw std::runtime_error(base_error);  // str contains invalid char(s)
}

inline auto from_chars(const char *first, const char *last, long double &value) noexcept {
  char **str_end = nullptr;
  long double buff{};
  if (*last == '\0') {  // strtold expects a null-terminated str
    buff = std::strtold(first, str_end);
  } else {
    const std::string strbuff{first, last};
    buff = std::strtold(strbuff.c_str(), str_end);
  }

  // NOLINTNEXTLINE(clang-analyzer-core.NullDereference)
  std::from_chars_result res{*str_end, std::errc()};
  if (res.ptr == first) {
    res.ec = std::errc::invalid_argument;
  } else if (buff == HUGE_VALL) {
    res.ec = std::errc::result_out_of_range;
  } else {
    value = buff;
  }

  return res;
}

template <typename N>
inline auto from_chars(const char *first, const char *last, N &value) noexcept {
  if constexpr (std::is_integral_v<N>) {
    return std::from_chars(first, last, value);
  } else {
    return fast_float::from_chars(first, last, value);
  }
}

template <typename N>
inline void parse_numeric_or_throw(std::string_view tok, N &field) {
  const auto *first = tok.data();
  const auto *last = first + tok.size();  // NOLINT
  auto [ptr, err] = from_chars(first, last, field);
  if (ptr != last || err != std::errc()) {
    throw_except_from_errc(tok, (std::numeric_limits<std::size_t>::max)(), field, ptr, err);
  }
}

template <typename N>
inline N parse_numeric_or_throw(std::string_view tok) {
  N field{};
  parse_numeric_or_throw(tok, field);
  return field;
}
}  // namespace hictk::internal
