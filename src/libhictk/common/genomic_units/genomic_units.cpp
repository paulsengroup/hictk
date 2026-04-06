// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/genomic_units.hpp"

#include <fmt/format.h>

#include <cctype>
#include <cstdint>
#include <stdexcept>
#include <string_view>

namespace hictk {

namespace internal {

void _throw_not_an_int_exception(std::string_view distance, double n) {
  throw std::runtime_error(
      fmt::format(FMT_STRING("Cannot convert \"{}\" ({} bp) to an integer number"), distance, n));
}

void _throw_int_out_of_bound_exception(double n, std::string_view type_name) {
  throw std::runtime_error(
      fmt::format(FMT_STRING("Cannot fit {:.0f} into a {} number"), n, type_name));
}

void _throw_genomic_unit_parse_error(std::string_view distance, const std::invalid_argument& e) {
  throw std::invalid_argument(fmt::format(
      FMT_STRING("failed to parse \"{}\" as genomic distance: {}"), distance, e.what()));
}

void _throw_genomic_unit_parse_error(std::string_view distance, const std::runtime_error& e) {
  throw std::runtime_error(fmt::format(FMT_STRING("failed to parse \"{}\" as genomic distance: {}"),
                                       distance, e.what()));
}

}  // namespace internal

std::uint32_t parse_genomic_unit(std::string_view unit) {
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

}  // namespace hictk
