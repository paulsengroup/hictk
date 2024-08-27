// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <string>
#include <toml++/toml.hpp>

#include "./common.hpp"
#include "./cool.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"

namespace hictk::tools {

[[nodiscard]] inline toml::table normalize_attribute_map(const cooler::SingleCellAttributes& map,
                                                         const std::string& uri) {
  toml::table attributes;

  if (!uri.empty()) {
    emplace_if_valid("uri", uri, attributes);
  }

  if (map.bin_size == 0) {
    assert(map.bin_type == BinTable::Type::variable);
    emplace_if_valid("bin-size", "variable", attributes);
  } else {
    assert(map.bin_type == BinTable::Type::fixed);
    emplace_if_valid("bin-size", map.bin_size, attributes);
  }
  emplace_if_valid("bin-type", map.bin_type == BinTable::Type::fixed ? "fixed" : "variable",
                   attributes);
  emplace_if_valid("format", map.format, attributes);
  emplace_if_valid("format-version", map.format_version, attributes);

  emplace_if_valid("creation-date", map.creation_date, attributes);
  emplace_if_valid("generated-by", map.generated_by, attributes);
  emplace_if_valid("assembly", map.assembly, attributes);
  emplace_if_valid("metadata", map.metadata, attributes);

  emplace_if_valid("format-url", map.format_url, attributes);
  emplace_if_valid("nbins", map.nbins, attributes);
  emplace_if_valid("ncells", map.ncells, attributes);
  emplace_if_valid("nchroms", map.nchroms, attributes);
  emplace_if_valid("storage-mode", map.storage_mode, attributes);

  return attributes;
}

[[nodiscard]] inline int print_scool_metadata(const std::filesystem::path& p,
                                              MetadataOutputFormat format, bool include_file_path,
                                              bool recursive) {
  const cooler::SingleCellFile sclr(p);
  auto attributes = normalize_attribute_map(sclr.attributes(), include_file_path ? p.string() : "");
  std::vector<std::pair<std::string, toml::table>> nested_attributes{};

  toml::array cells;
  for (const auto& cell : sclr.cells()) {
    cells.push_back(cell);
  }
  emplace_if_valid("cells", cells, attributes);

  if (recursive) {
    for (const auto& cell_id : sclr.cells()) {
      nested_attributes.emplace_back(cell_id,
                                     normalize_attribute_map(sclr.open(cell_id).attributes(), ""));
    }
  }

  print_attributes(attributes, nested_attributes, format);
  return 0;
}

}  // namespace hictk::tools
