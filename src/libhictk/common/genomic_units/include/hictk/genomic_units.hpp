// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <type_traits>

#include "hictk/common.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/type_pretty_printer.hpp"

namespace hictk {

namespace internal {

[[noreturn]] void _throw_not_an_int_exception(std::string_view distance, double n);
[[noreturn]] void _throw_int_out_of_bound_exception(double n, std::string_view type_name);
[[noreturn]] void _throw_genomic_unit_parse_error(std::string_view distance,
                                                  const std::invalid_argument& e);
[[noreturn]] void _throw_genomic_unit_parse_error(std::string_view distance,
                                                  const std::runtime_error& e);

}  // namespace internal

[[nodiscard]] std::uint32_t parse_genomic_unit(std::string_view unit);

template <typename T = std::uint64_t>
[[nodiscard]] T parse_genomic_distance(std::string_view distance) {
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
        internal::_throw_not_an_int_exception(distance, n);
      }
      if (n < static_cast<double>(std::numeric_limits<T>::lowest()) ||
          n > static_cast<double>(std::numeric_limits<T>::max())) {
        internal::_throw_int_out_of_bound_exception(n, internal::type_name<T>());
      }
    }
    return conditional_static_cast<T>(n);
  } catch (const std::invalid_argument& e) {
    internal::_throw_genomic_unit_parse_error(distance, e);
  } catch (const std::runtime_error& e) {
    internal::_throw_genomic_unit_parse_error(distance, e);
  }
}

}  // namespace hictk
