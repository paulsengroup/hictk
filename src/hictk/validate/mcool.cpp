// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include "./validate.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/tools/file_attributes_formatting.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {

static void update_status_table(const cooler::utils::ValidationStatusMultiresCooler& status,
                                toml::table& buff) {
  buff.insert("is_hdf5", status.is_hdf5);
  buff.insert("unable_to_open_file", status.unable_to_open_file);
  buff.insert("file_was_properly_closed", status.file_was_properly_closed);
  buff.insert("missing_or_invalid_format_attr", status.missing_or_invalid_format_attr);
  buff.insert("missing_or_invalid_bin_type_attr", status.missing_or_invalid_bin_type_attr);
  buff.insert("missing_groups", io::toml::to_array(status.missing_groups));
  buff.insert("is_valid_mcool", status.is_multires_file);

  assert(status.valid_resolutions.empty());
  assert(status.invalid_resolutions.empty());
}

[[nodiscard]] static std::optional<cooler::MultiResFile> open_mcool_noexcept(
    std::string_view uri) noexcept {
  try {
    return cooler::MultiResFile{std::filesystem::path(uri)};
  } catch (const std::exception& e) {
    SPDLOG_DEBUG(FMT_STRING("failed to open file \"{}\": {}"), uri, e.what());
    return {};
  } catch (...) {
    SPDLOG_DEBUG(FMT_STRING("failed to open file \"{}\": unknown error"), uri);
    return {};
  }
}

// NOLINTNEXTLINE(bugprone-exception-escape)
[[nodiscard]] static std::string get_cooler_uri_noexcept(const cooler::MultiResFile& mclr,
                                                         std::uint32_t resolution) noexcept {
  try {
    return std::string{mclr.open(resolution).uri()};
  } catch ([[maybe_unused]] const std::exception& e) {
    SPDLOG_DEBUG(FMT_STRING("failed to open Cooler at resolution {} from file \"{}\": {}"),
                 resolution, mclr.path(), e.what());
  } catch (...) {
    SPDLOG_DEBUG(
        FMT_STRING("failed to open Cooler at resolution {} from file \"{}\": unknown error"),
        resolution, mclr.path());
  }

  return fmt::format(FMT_STRING("{}::/resolutions/{}"), mclr.path(), resolution);
}

std::pair<int, toml::table> validate_mcool(std::string_view path, bool validate_index,
                                           bool validate_pixels, bool exhaustive) {
  int return_code = 0;
  toml::table global_status;
  const auto validation_status = cooler::utils::is_multires_file(path, false);
  update_status_table(validation_status, global_status);

  if (!validation_status.is_multires_file) {
    return std::make_pair(1, global_status);
  }

  const auto mclr = open_mcool_noexcept(path);
  if (!mclr) {
    global_status.insert_or_assign("is_valid_mcool", false);
    return std::make_pair(1, global_status);
  }

  bool early_return = false;
  std::for_each(mclr->resolutions().rbegin(), mclr->resolutions().rend(), [&](const auto res) {
    if (early_return) {
      return;
    }
    const auto [_, status] =
        validate_cooler(get_cooler_uri_noexcept(*mclr, res), validate_index, validate_pixels);
    const auto is_cooler = **status.template get_as<bool>("is_valid_cooler");
    global_status.insert(fmt::to_string(res), status);
    if (!is_cooler) {
      return_code = 1;
      if (!exhaustive) {
        early_return = true;
      }
    }
  });

  if (return_code != 0) {
    global_status.insert_or_assign("is_valid_mcool", false);
  }

  return std::make_pair(return_code, global_status);
}

}  // namespace hictk::tools
