// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <filesystem>
#include <string>

#include "./metadata.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {

toml::table normalize_attribute_map(const cooler::Attributes& map, const std::string& uri) {
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
  emplace_if_valid("storage-mode", map.storage_mode, attributes);

  emplace_if_valid("creation-date", map.creation_date, attributes);
  emplace_if_valid("generated-by", map.generated_by, attributes);
  emplace_if_valid("assembly", map.assembly, attributes);
  emplace_if_valid("metadata", map.metadata, attributes);

  emplace_if_valid("format-url", map.format_url, attributes);
  emplace_if_valid("nbins", map.nbins, attributes);
  emplace_if_valid("nchroms", map.nchroms, attributes);
  emplace_if_valid("nnz", map.nnz, attributes);
  emplace_if_valid("sum", map.sum, attributes);
  emplace_if_valid("cis", map.cis, attributes);

  return attributes;
}

int print_cool_metadata(const std::filesystem::path& p, MetadataOutputFormat format,
                        bool include_file_path) {
  const auto attributes = normalize_attribute_map(cooler::File(p.string()).attributes(),
                                                  include_file_path ? p.string() : "");
  print_attributes(attributes, {}, format);

  return 0;
}

}  // namespace hictk::tools
