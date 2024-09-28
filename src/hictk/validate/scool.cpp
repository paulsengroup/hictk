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
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/tools/file_attributes_formatting.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {

static void update_status_table(const cooler::utils::ValidationStatusScool& status,
                                toml::table& buff) {
  buff.insert("is_hdf5", status.is_hdf5);
  buff.insert("unable_to_open_file", status.unable_to_open_file);
  buff.insert("file_was_properly_closed", status.file_was_properly_closed);
  buff.insert("missing_or_invalid_format_attr", status.missing_or_invalid_format_attr);
  buff.insert("missing_or_invalid_bin_type_attr", status.missing_or_invalid_bin_type_attr);
  buff.insert("missing_groups", io::toml::to_array(status.missing_groups));
  buff.insert("unexpected_number_of_cells", status.unexpected_number_of_cells);
  buff.insert("is_valid_scool", status.is_scool_file);

  assert(status.valid_cells.empty());
  assert(status.invalid_cells.empty());
}

[[nodiscard]] static std::optional<cooler::SingleCellFile> open_scool_noexcept(
    std::string_view uri) noexcept {
  try {
    return cooler::SingleCellFile{std::filesystem::path(uri)};
  } catch (const std::exception& e) {
    SPDLOG_DEBUG(FMT_STRING("failed to open file \"{}\": {}"), uri, e.what());
    return {};
  } catch (...) {
    SPDLOG_DEBUG(FMT_STRING("failed to open file \"{}\": unknown error"), uri);
    return {};
  }
}

std::pair<int, toml::table> validate_scool(std::string_view path, bool validate_index,
                                           bool exhaustive) {
  int return_code = 0;
  toml::table global_status;

  update_status_table(cooler::utils::is_scool_file(path, false), global_status);
  const auto sclr = open_scool_noexcept(path);
  if (!sclr) {
    global_status.insert_or_assign("is_valid_scool", false);
    return std::make_pair(1, global_status);
  }

  for (const auto& cell : sclr->cells()) {
    const auto [_, status] = validate_cooler(sclr->open(cell).uri(), validate_index);
    const auto is_cooler = **status.get_as<bool>("is_valid_cooler");
    global_status.insert(cell, status);

    if (!is_cooler) {
      return_code = 1;
      if (!exhaustive) {
        global_status.insert_or_assign("is_valid_scool", false);
        return std::make_pair(1, global_status);
      }
    }
  }

  if (return_code != 0) {
    global_status.insert_or_assign("is_valid_scool", false);
  }

  return std::make_pair(return_code, global_status);
}

}  // namespace hictk::tools
