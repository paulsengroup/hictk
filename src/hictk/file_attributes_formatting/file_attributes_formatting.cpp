// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/file_attributes_formatting.hpp"

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <limits>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/tools/toml.hpp"
#include "hictk/version.hpp"

namespace hictk::tools::io {

namespace json {

[[nodiscard]] nlohmann::json reformat_nulls(nlohmann::json attributes) {
  std::vector<std::string> null_fields{};
  for (const auto& field : attributes.items()) {
    if (field.value() == "null") {
      null_fields.emplace_back(field.key());
    }
  }

  for (const auto& k : null_fields) {
    attributes[k] = nullptr;
  }

  return attributes;
}

}  // namespace json

namespace toml {

[[nodiscard]] nlohmann::json toml_to_json(const ::toml::table& t) {
  std::stringstream buff;

  // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
  buff << ::toml::json_formatter(t);

  auto j = json::reformat_nulls(nlohmann::json::parse(buff.str()));
  if (const auto metadata = t.find("metadata");
      metadata != t.end() && metadata->second.is_string()) {
    try {
      j["metadata"] =
          json::reformat_nulls(nlohmann::json::parse(metadata->second.ref<std::string>()));
      // NOLINTNEXTLINE
    } catch (...) {
    }
  }

  return j;
}

[[nodiscard]] std::string format_to_json(
    const ::toml::table& attributes,
    const std::vector<std::pair<std::string, ::toml::table>>& nested_attributes) {
  // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
  auto attributes_json = toml_to_json(attributes);

  for (const auto& [key, table] : nested_attributes) {
    attributes_json[key] = toml_to_json(table);
  }

  return attributes_json.dump(4);
}

[[nodiscard]] std::string sanitize_toml_section_title(std::string s) {
  if (s.find('.') == std::string::npos) {
    return s;
  }

  // Escape '
  std::size_t start_pos = 0;
  while ((start_pos = s.find('\'', start_pos)) != std::string::npos) {
    s.replace(start_pos, 1, "\\'");
    start_pos += 2;
  }

  s.insert(s.begin(), 1, '\'');
  s.insert(s.end(), 1, '\'');

  return s;
}

[[nodiscard]] std::string format_to_toml(
    const ::toml::table& attributes,
    const std::vector<std::pair<std::string, ::toml::table>>& nested_attributes) {
  std::stringstream ss;
  ss << fmt::format(FMT_STRING("# Metadata generated by {}\n"), hictk::config::version::str_long())
     << attributes << '\n';

  for (const auto& [title, table] : nested_attributes) {
    ss << fmt::format(FMT_STRING("\n[{}]\n"), sanitize_toml_section_title(title)) << table << '\n';
  }

  return ss.str();
}

[[nodiscard]] std::string format_to_yaml(
    const ::toml::table& attributes,
    const std::vector<std::pair<std::string, ::toml::table>>& nested_attributes) {
  // NOLINTBEGIN(clang-analyzer-optin.core.EnumCastOutOfRange)
  if (nested_attributes.empty()) {
    std::stringstream ss;
    ss << fmt::format(FMT_STRING("--- # Metadata generated by {}\n"),
                      hictk::config::version::str_long())

       << ::toml::yaml_formatter(attributes) << '\n';
    return ss.str();
  }
  // NOLINTEND(clang-analyzer-optin.core.EnumCastOutOfRange)

  return format_to_yaml(::toml::parse(format_to_toml(attributes, nested_attributes)), {});
}
}  // namespace toml
}  // namespace hictk::tools::io
