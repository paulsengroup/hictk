// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
// clang-format off
#include "hictk/suppress_warnings.hpp"
// clang-format on
#include <fmt/format.h>

#include <cstdint>
DISABLE_WARNING_PUSH
DISABLE_WARNING_NULL_DEREF
#include <highfive/H5File.hpp>
DISABLE_WARNING_POP
#include <highfive/H5Group.hpp>
#include <string_view>
#include <utility>

namespace hictk::cooler::utils {

namespace internal {
struct ValidationStatusBase {
  bool is_hdf5{false};
  bool file_was_properly_closed{false};

  bool missing_or_invalid_format_attr{true};
  bool missing_or_invalid_bin_type_attr{true};

  std::string uri{};
  std::vector<std::string> missing_groups{};
};
}  // namespace internal

struct ValidationStatusCooler : public internal::ValidationStatusBase {
  bool is_cooler{false};

  constexpr explicit operator bool() const noexcept;
};

struct ValidationStatusMultiresCooler : public internal::ValidationStatusBase {
  bool is_multires_file{false};

  std::vector<ValidationStatusCooler> invalid_resolutions{};

  constexpr explicit operator bool() const noexcept;
};

struct ValidationStatusScool : public internal::ValidationStatusBase {
  bool is_scool_file{false};

  bool unexpected_number_of_cells{true};
  std::vector<ValidationStatusCooler> invalid_cells{};

  constexpr explicit operator bool() const noexcept;
};

[[nodiscard]] ValidationStatusCooler is_cooler(std::string_view uri);
[[nodiscard]] ValidationStatusCooler is_cooler(const HighFive::File& fp,
                                               std::string_view root_path = "/");
[[nodiscard]] ValidationStatusCooler is_cooler(const HighFive::Group& root_group);

[[nodiscard]] ValidationStatusMultiresCooler is_multires_file(std::string_view uri,
                                                              bool validate_resolutions = true,
                                                              std::int64_t min_version = 1);
[[nodiscard]] ValidationStatusMultiresCooler is_multires_file(const HighFive::File& fp,
                                                              bool validate_resolutions = true,
                                                              std::int64_t min_version = 1);

[[nodiscard]] ValidationStatusScool is_scool_file(std::string_view uri, bool validate_cells = true);
[[nodiscard]] ValidationStatusScool is_scool_file(const HighFive::File& fp,
                                                  bool validate_cells = true);

[[nodiscard]] bool index_is_valid(std::string_view uri, bool verbose = false);
}  // namespace hictk::cooler::utils

namespace fmt {
template <>
struct formatter<hictk::cooler::utils::ValidationStatusCooler> {
  constexpr auto parse(format_parse_context& ctx) const -> format_parse_context::iterator;

  inline auto format(const hictk::cooler::utils::ValidationStatusCooler& s,
                     format_context& ctx) const -> format_context::iterator;
};

template <>
struct formatter<hictk::cooler::utils::ValidationStatusMultiresCooler> {
  constexpr auto parse(format_parse_context& ctx) const -> format_parse_context::iterator;

  inline auto format(const hictk::cooler::utils::ValidationStatusMultiresCooler& s,
                     format_context& ctx) const -> format_context::iterator;
};

template <>
struct formatter<hictk::cooler::utils::ValidationStatusScool> {
  constexpr auto parse(format_parse_context& ctx) const -> format_parse_context::iterator;

  inline auto format(const hictk::cooler::utils::ValidationStatusScool& s,
                     format_context& ctx) const -> format_context::iterator;
};
}  // namespace fmt

#include "../../../validation_impl.hpp"
