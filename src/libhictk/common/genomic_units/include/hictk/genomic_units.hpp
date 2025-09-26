// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <type_traits>

#include "hictk/common.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/type_pretty_printer.hpp"

namespace hictk {

[[nodiscard]] inline std::uint32_t parse_genomic_unit(std::string_view unit) {
  auto handle_invalid_unit = [&]() {
    return std::invalid_argument(
        fmt::format(FMT_STRING("Unrecognized unit \"{}\": valid units are k[bp], m[bp], and "
                               "g[bp] (case-insensitive, e.g. k, KB, or KBP)."),
                    unit));
  };

  if (unit.empty()) {
    throw std::invalid_argument("unit is empty");
  }

  if (unit.size() == 2 && std::tolower(unit.front()) == 'b' && std::tolower(unit.back()) == 'p') {
    return 1;
  }

  if (unit.size() > 3) {
    throw handle_invalid_unit();
  }
  if (unit.size() > 1 && std::tolower(unit[1]) != 'b') {
    throw handle_invalid_unit();
  }
  if (unit.size() > 2 && std::tolower(unit[2]) != 'p') {
    throw handle_invalid_unit();
  }

  // NOLINTBEGIN(*-avoid-magic-numbers)
  switch (std::tolower(unit.front())) {
    case 'k':
      return 1'000;
    case 'm':
      return 1'000'000;
    case 'g':
      return 1'000'000'000;
    default:
      throw handle_invalid_unit();
  }
  // NOLINTEND(*-avoid-magic-numbers)
}

template <typename T = std::uint64_t>
[[nodiscard]] inline T parse_genomic_distance(std::string_view distance) {
  static_assert(std::is_arithmetic_v<T>);

  try {
    if (distance.empty()) {
      throw std::invalid_argument("distance is empty");
    }

    auto match = std::find_if(distance.begin(), distance.end(),
                              [](char c) { return !std::isdigit(c) && c != '.'; });

    if (match == distance.begin() || distance.front() == '.') {
      throw std::invalid_argument("distance does not start with a digit");
    }

    if (match == distance.end()) {
      return internal::parse_numeric_or_throw<T>(distance);
    }

    auto size = static_cast<std::size_t>(std::distance(distance.begin(), match));
    const auto cfx = distance.substr(0, size);

    match = std::find_if(match, distance.end(), [](char c) { return !std::isspace(c); });

    if (match == distance.end()) {
      throw std::invalid_argument("distance has trailing whitespaces");
    }

    size = static_cast<std::size_t>(std::distance(distance.begin(), match));
    const auto unit = distance.substr(size);

    const auto numeric_cfx = internal::parse_numeric_or_throw<double>(cfx);
    const auto multiplier = static_cast<double>(parse_genomic_unit(unit));
    const auto n = numeric_cfx * multiplier;
    if constexpr (!std::is_floating_point_v<T>) {
      if (n != std::floor(n)) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("Cannot convert \"{}\" ({} bp) to an integer number"), distance, n));
      }
      if (n < static_cast<double>(std::numeric_limits<T>::lowest()) ||
          n > static_cast<double>(std::numeric_limits<T>::max())) {
        throw std::runtime_error(fmt::format(FMT_STRING("Cannot fit {:.0f} into a {} number"), n,
                                             internal::type_name<T>()));
      }
    }
    return conditional_static_cast<T>(n);
  } catch (const std::invalid_argument& e) {
    throw std::invalid_argument(fmt::format(
        FMT_STRING("failed to parse \"{}\" as genomic distance: {}"), distance, e.what()));
  } catch (const std::runtime_error& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("failed to parse \"{}\" as genomic distance: {}"), distance, e.what()));
  }
}

}  // namespace hictk
