// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <nlohmann/json.hpp>
#include <optional>
#include <sstream>
#include <string>
#include <string_view>
#include <toml++/toml.hpp>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/version.hpp"

namespace hictk::tools {
enum class MetadataOutputFormat : std::uint8_t { json, toml, yaml };

[[nodiscard]] inline MetadataOutputFormat parse_output_format(std::string_view format) {
  if (format == "json") {
    return MetadataOutputFormat::json;
  }
  if (format == "toml") {
    return MetadataOutputFormat::toml;
  }
  assert(format == "yaml");
  return MetadataOutputFormat::yaml;
}

inline void emplace_if_valid(std::string_view key, const std::string& value, toml::table& buff) {
  if (!key.empty()) {
    buff.insert(key, value);
  }
}

inline void emplace_if_valid(std::string_view key, const toml::array& values, toml::table& buff) {
  if (!key.empty()) {
    buff.insert(key, values);
  }
}

template <typename T, typename std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline void emplace_if_valid(std::string_view key, const T& value, toml::table& buff) {
  if (key.empty()) {
    return;
  }

  if (value <= std::numeric_limits<std::int64_t>::max()) {
    buff.insert(key, static_cast<std::int64_t>(value));
  } else {
    emplace_if_valid(key, fmt::to_string(value), buff);
  }
}

template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>>* = nullptr>
inline void emplace_if_valid(std::string_view key, const T& value, toml::table& buff) {
  if (!key.empty()) {
    buff.insert(key, conditional_static_cast<double>(value));
  }
}

template <typename... T>
inline void emplace_if_valid(std::string_view key, const std::variant<T...>& value,
                             toml::table& buff) {
  if (!key.empty()) {
    std::visit([&](const auto& value_) { emplace_if_valid(key, value_, buff); }, value);
  }
}

template <typename T>
inline void emplace_if_valid(std::string_view key, const std::optional<T>& value,
                             toml::table& buff) {
  if (!key.empty() && !!value) {
    emplace_if_valid(key, *value, buff);
  }
}

[[nodiscard]] inline nlohmann::json reformat_nulls(nlohmann::json attributes) {
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

[[nodiscard]] inline nlohmann::json toml_to_json(const toml::table& t) {
  std::stringstream buff;
  buff << toml::json_formatter(t);

  auto j = reformat_nulls(nlohmann::json::parse(buff.str()));
  if (const auto metadata = t.find("metadata");
      metadata != t.end() && metadata->second.is_string()) {
    try {
      j["metadata"] = reformat_nulls(nlohmann::json::parse(metadata->second.ref<std::string>()));
      // NOLINTNEXTLINE
    } catch (...) {
    }
  }

  return j;
}

[[nodiscard]] inline std::string format_to_json(
    const toml::table& attributes,
    const std::vector<std::pair<std::string, toml::table>>& nested_attributes) {
  auto attributes_json = toml_to_json(attributes);

  for (const auto& [key, table] : nested_attributes) {
    attributes_json[key] = toml_to_json(table);
  }

  return attributes_json.dump(4);
}

[[nodiscard]] inline std::string sanitize_toml_section_title(std::string s) {
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

[[nodiscard]] inline std::string format_to_toml(
    const toml::table& attributes,
    const std::vector<std::pair<std::string, toml::table>>& nested_attributes) {
  std::stringstream ss;
  ss << fmt::format(FMT_STRING("# Metadata generated by {}\n"), hictk::config::version::str_long())
     << attributes << '\n';

  for (const auto& [title, table] : nested_attributes) {
    ss << fmt::format(FMT_STRING("\n[{}]\n"), sanitize_toml_section_title(title)) << table << '\n';
  }

  return ss.str();
}

[[nodiscard]] inline std::string format_to_yaml(
    const toml::table& attributes,
    const std::vector<std::pair<std::string, toml::table>>& nested_attributes) {
  if (nested_attributes.empty()) {
    std::stringstream ss;
    ss << fmt::format(FMT_STRING("--- # Metadata generated by {}\n"),
                      hictk::config::version::str_long())
       << toml::yaml_formatter(attributes) << '\n';
    return ss.str();
  }

  const auto metadata_toml = toml::parse(format_to_toml(attributes, nested_attributes));

  return format_to_yaml(metadata_toml, {});
}

inline void print_attributes(
    const toml::table& top_lvl_attributes,
    const std::vector<std::pair<std::string, toml::table>>& nested_attributes,
    MetadataOutputFormat format) {
  std::string buff;
  switch (format) {
    case MetadataOutputFormat::json: {
      buff = format_to_json(top_lvl_attributes, nested_attributes);
      break;
    }
    case MetadataOutputFormat::toml: {
      buff = format_to_toml(top_lvl_attributes, nested_attributes);
      break;
    }
    case MetadataOutputFormat::yaml: {
      buff = format_to_yaml(top_lvl_attributes, nested_attributes);
      break;
    }
  }

  fmt::print(FMT_STRING("{}\n"), buff);
}

}  // namespace hictk::tools
