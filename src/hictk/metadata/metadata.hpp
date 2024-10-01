// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <filesystem>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {
enum class MetadataOutputFormat : std::uint8_t { json, toml, yaml };

[[nodiscard]] toml::table normalize_attribute_map(const cooler::Attributes& map,
                                                  const std::string& uri);

[[nodiscard]] MetadataOutputFormat parse_output_format(std::string_view format);

void emplace_if_valid(std::string_view key, const std::string& value, toml::table& buff);

void emplace_if_valid(std::string_view key, const toml::array& values, toml::table& buff);

template <typename T, typename std::enable_if_t<std::is_integral_v<T>>* = nullptr>
inline void emplace_if_valid(std::string_view key, const T& value, toml::table& buff) {
  if (key.empty()) {
    return;
  }
  HICTK_DISABLE_WARNING_PUSH
  HICTK_DISABLE_WARNING_BOOL_COMPARE
  if (value <= std::numeric_limits<std::int64_t>::max()) {
    buff.insert(key, static_cast<std::int64_t>(value));
  } else {
    emplace_if_valid(key, fmt::to_string(value), buff);
  }
  HICTK_DISABLE_WARNING_POP
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

[[nodiscard]] int print_cool_metadata(const std::filesystem::path& p, MetadataOutputFormat format,
                                      bool include_file_path);

[[nodiscard]] int print_hic_metadata(const std::filesystem::path& p, MetadataOutputFormat format,
                                     bool include_file_path, bool recursive);

[[nodiscard]] int print_mcool_metadata(const std::filesystem::path& p, MetadataOutputFormat format,
                                       bool include_file_path, bool recursive);

[[nodiscard]] int print_scool_metadata(const std::filesystem::path& p, MetadataOutputFormat format,
                                       bool include_file_path, bool recursive);

void print_attributes(const toml::table& top_lvl_attributes,
                      const std::vector<std::pair<std::string, toml::table>>& nested_attributes,
                      MetadataOutputFormat format);

}  // namespace hictk::tools
