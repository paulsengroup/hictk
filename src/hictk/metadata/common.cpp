// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <cassert>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "./metadata.hpp"
#include "hictk/tools/file_attributes_formatting.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {

[[nodiscard]] MetadataOutputFormat parse_output_format(std::string_view format) {
  if (format == "json") {
    return MetadataOutputFormat::json;
  }
  if (format == "toml") {
    return MetadataOutputFormat::toml;
  }
  assert(format == "yaml");
  return MetadataOutputFormat::yaml;
}

void emplace_if_valid(std::string_view key, const std::string& value, toml::table& buff) {
  if (!key.empty()) {
    buff.insert(key, value);
  }
}

void emplace_if_valid(std::string_view key, const toml::array& values, toml::table& buff) {
  if (!key.empty()) {
    buff.insert(key, values);
  }
}

void print_attributes(const toml::table& top_lvl_attributes,
                      const std::vector<std::pair<std::string, toml::table>>& nested_attributes,
                      MetadataOutputFormat format) {
  std::string buff;
  switch (format) {
    case MetadataOutputFormat::json: {
      buff = io::toml::format_to_json(top_lvl_attributes, nested_attributes);
      break;
    }
    case MetadataOutputFormat::toml: {
      buff = io::toml::format_to_toml(top_lvl_attributes, nested_attributes);
      break;
    }
    case MetadataOutputFormat::yaml: {
      buff = io::toml::format_to_yaml(top_lvl_attributes, nested_attributes);
      break;
    }
  }

  fmt::print(FMT_STRING("{}\n"), buff);
}

}  // namespace hictk::tools
