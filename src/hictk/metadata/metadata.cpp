// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <nlohmann/json.hpp>
#include <string>
#include <string_view>
#include <toml.hpp>
#include <utility>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/version.hpp"

namespace hictk::tools {

enum class MetadataOutputFormat : std::uint8_t { json, toml, yaml };

using AttributeValue = std::variant<std::int64_t, double, bool, std::string>;

[[nodiscard]] static AttributeValue try_parse_str(const std::string& value) {
  try {
    // NOLINTNEXTLINE
    return {internal::parse_numeric_or_throw<std::int64_t>(value)};
  } catch (...) {  // NOLINT
  }

  try {
    // NOLINTNEXTLINE
    return {internal::parse_numeric_or_throw<double>(value)};
  } catch (...) {  // NOLINT
  }

  if (value == "true" || value == "True") {
    return {true};
  }
  if (value == "false" || value == "False") {
    return {false};
  }

  if (value == "NULL" || value == "Null" || value == "null" || value == "None") {
    return {"null"};
  }

  return {value};
}

template <typename... T>
[[nodiscard]] static AttributeValue try_parse_str(const std::variant<T...>& var) {
  return std::visit([&](const auto& value) { return try_parse_str(value); }, var);
}

template <typename T>
[[nodiscard]] static AttributeValue try_parse_str(const std::optional<T>& opt) {
  if (opt.has_value()) {
    return try_parse_str(opt.value());
  }
  return {"null"};
}

[[nodiscard]] static MetadataOutputFormat parse_output_format(std::string_view format) {
  if (format == "json") {
    return MetadataOutputFormat::json;
  }
  if (format == "toml") {
    return MetadataOutputFormat::toml;
  }
  assert(format == "yaml");
  return MetadataOutputFormat::yaml;
}

static void emplace_if_valid(std::string_view key, const std::string& value, toml::table& buff) {
  if (!key.empty()) {
    buff.insert(key, value);
  }
}

template <typename T, typename std::enable_if_t<std::is_integral_v<T>>* = nullptr>
static void emplace_if_valid(std::string_view key, const T& value, toml::table& buff) {
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
static void emplace_if_valid(std::string_view key, const T& value, toml::table& buff) {
  if (!key.empty()) {
    buff.insert(key, conditional_static_cast<double>(value));
  }
}

template <typename... T>
static void emplace_if_valid(std::string_view key, const std::variant<T...>& value,
                             toml::table& buff) {
  if (!key.empty()) {
    std::visit([&](const auto& value_) { emplace_if_valid(key, value_, buff); }, value);
  }
}

template <typename T>
static void emplace_if_valid(std::string_view key, const std::optional<T>& value,
                             toml::table& buff) {
  if (!key.empty() && !!value) {
    emplace_if_valid(key, *value, buff);
  }
}

[[nodiscard]] static toml::table normalize_attribute_map(
    const phmap::flat_hash_map<std::string, std::string>& map, const std::string& uri) {
  toml::table attributes;

  if (!uri.empty()) {
    emplace_if_valid("uri", uri, attributes);
  }

  for (const auto& [k, v] : map) {
    std::visit([&, key = k](const auto& value) { emplace_if_valid(key, value, attributes); },
               try_parse_str(v));
  }

  return attributes;
}

[[nodiscard]] static toml::table normalize_attribute_map(const cooler::MultiResAttributes& map,
                                                         const std::string& uri) {
  toml::table attributes;

  if (!uri.empty()) {
    emplace_if_valid("uri", uri, attributes);
  }

  emplace_if_valid("bin-type", map.bin_type, attributes);
  emplace_if_valid("format", map.format, attributes);
  emplace_if_valid("format-version", map.format_version, attributes);

  return attributes;
}

[[nodiscard]] static toml::table normalize_attribute_map(const cooler::SingleCellAttributes& map,
                                                         const std::string& uri) {
  toml::table attributes;

  if (!uri.empty()) {
    emplace_if_valid("uri", uri, attributes);
  }

  emplace_if_valid("bin-size", map.bin_size, attributes);
  emplace_if_valid("bin-type", map.bin_type, attributes);
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

[[nodiscard]] static toml::table normalize_attribute_map(const cooler::Attributes& map,
                                                         const std::string& uri) {
  toml::table attributes;

  if (!uri.empty()) {
    emplace_if_valid("uri", uri, attributes);
  }

  emplace_if_valid("bin-size", map.bin_size, attributes);
  emplace_if_valid("bin-type", map.bin_type, attributes);
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

[[nodiscard]] static nlohmann::json reformat_nulls(nlohmann::json attributes) {
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

static void format_to_json(const toml::table& attributes) {
  if (!attributes.contains("metadata") || !attributes["metadata"].is_string()) {
    std::cout << toml::json_formatter(attributes) << '\n';
    return;
  }

  try {
    // Try to pretty-print metadata attribute
    auto new_attributes = attributes;
    new_attributes.erase(new_attributes.find("metadata"));
    std::stringstream buff;
    buff << toml::json_formatter(new_attributes);

    auto attributes_json = nlohmann::json::parse(buff.str());
    attributes_json["metadata"] =
        reformat_nulls(nlohmann::json::parse(attributes["metadata"].ref<std::string>()));

    fmt::print(FMT_STRING("{}\n"), attributes_json.dump(4));

  } catch (...) {
    std::cout << toml::json_formatter(attributes) << '\n';
  }
}

static void format_to_toml(const toml::table& attributes) {
  std::cout << fmt::format(FMT_STRING("# Metadata generated by {}\n"),
                           hictk::config::version::str_long())
            << attributes << '\n';
}

static void format_to_yaml(const toml::table& attributes) {
  std::cout << fmt::format(FMT_STRING("--- # Metadata generated by {}\n"),
                           hictk::config::version::str_long())
            << toml::yaml_formatter(attributes) << '\n';
}

static void print_attributes(const toml::table& attributes, MetadataOutputFormat format) {
  switch (format) {
    case MetadataOutputFormat::json:
      return format_to_json(attributes);  // NOLINT
    case MetadataOutputFormat::toml:
      return format_to_toml(attributes);  // NOLINT
    case MetadataOutputFormat::yaml:
      return format_to_yaml(attributes);  // NOLINT
  }
}

[[nodiscard]] static int print_hic_metadata(const std::filesystem::path& p,
                                            MetadataOutputFormat format, bool include_file_path) {
  const auto resolution = hic::utils::list_resolutions(p).front();
  const auto attributes = normalize_attribute_map(hic::File(p, resolution).attributes(),
                                                  include_file_path ? p.string() : "");
  print_attributes(attributes, format);

  return 0;
}

[[nodiscard]] static int print_mcool_metadata(const std::filesystem::path& p,
                                              MetadataOutputFormat format, bool include_file_path) {
  const auto attributes = normalize_attribute_map(cooler::MultiResFile(p).attributes(),
                                                  include_file_path ? p.string() : "");
  print_attributes(attributes, format);

  return 0;
}

[[nodiscard]] static int print_scool_metadata(const std::filesystem::path& p,
                                              MetadataOutputFormat format, bool include_file_path) {
  const auto attributes = normalize_attribute_map(cooler::SingleCellFile(p).attributes(),
                                                  include_file_path ? p.string() : "");
  print_attributes(attributes, format);

  return 0;
}

[[nodiscard]] static int print_cool_metadata(const std::filesystem::path& p,
                                             MetadataOutputFormat format, bool include_file_path) {
  const auto attributes = normalize_attribute_map(cooler::File(p.string()).attributes(),
                                                  include_file_path ? p.string() : "");
  print_attributes(attributes, format);

  return 0;
}

int metadata_subcmd(const MetadataConfig& c) {
  const auto output_format = parse_output_format(c.output_format);
  if (c.input_format == "hic") {
    return print_hic_metadata(c.uri, output_format, c.include_file_path);
  }
  if (c.input_format == "mcool") {
    return print_mcool_metadata(c.uri, output_format, c.include_file_path);
  }
  if (c.input_format == "scool") {
    return print_scool_metadata(c.uri, output_format, c.include_file_path);
  }
  assert(c.input_format == "cool");
  return print_cool_metadata(c.uri, output_format, c.include_file_path);
}

}  // namespace hictk::tools
