// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "./metadata.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {

[[nodiscard]] static toml::table normalize_attribute_map(const cooler::MultiResAttributes& map,
                                                         const std::string& uri) {
  toml::table attributes;

  if (!uri.empty()) {
    emplace_if_valid("uri", uri, attributes);
  }

  emplace_if_valid("bin-type", map.bin_type == BinTable::Type::fixed ? "fixed" : "variable",
                   attributes);
  emplace_if_valid("format", map.format, attributes);
  emplace_if_valid("format-version", map.format_version, attributes);

  return attributes;
}

int print_mcool_metadata(const std::filesystem::path& p, MetadataOutputFormat format,
                         bool include_file_path, bool recursive) {
  const cooler::MultiResFile mclr(p);
  auto attributes = normalize_attribute_map(mclr.attributes(), include_file_path ? p.string() : "");
  std::vector<std::pair<std::string, toml::table>> nested_attributes{};

  toml::array resolutions;
  for (const auto& resolution : mclr.resolutions()) {
    if (resolution == 0) {
      resolutions.push_back("variable");
    } else {
      resolutions.push_back(static_cast<std::int64_t>(resolution));
    }
  }
  emplace_if_valid("resolutions", resolutions, attributes);

  if (recursive) {
    for (const auto& resolution : mclr.resolutions()) {
      nested_attributes.emplace_back(
          fmt::to_string(resolution),
          normalize_attribute_map(mclr.open(resolution).attributes(), ""));
    }
  }

  print_attributes(attributes, nested_attributes, format);
  return 0;
}

}  // namespace hictk::tools
