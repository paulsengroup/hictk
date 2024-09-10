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
#include <string>
#include <string_view>
#include <toml++/toml.hpp>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/tools/common.hpp"

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
  DISABLE_WARNING_PUSH
  DISABLE_WARNING_BOOL_COMPARE
  if (value <= std::numeric_limits<std::int64_t>::max()) {
    buff.insert(key, static_cast<std::int64_t>(value));
  } else {
    emplace_if_valid(key, fmt::to_string(value), buff);
  }
  DISABLE_WARNING_POP
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

inline void print_attributes(
    const toml::table& top_lvl_attributes,
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
