// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic.hpp"

#include <fmt/format.h>

#include <cstdint>
#include <filesystem>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include "./metadata.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/tools/toml.hpp"

namespace hictk::tools {
using AttributeValue = std::variant<std::int64_t, double, bool, std::string>;

[[nodiscard]] static AttributeValue try_parse_str(const std::string& value) {
  try {
    return {internal::parse_numeric_or_throw<std::int64_t>(value)};
    // NOLINTNEXTLINE
  } catch (...) {
  }

  try {
    return {internal::parse_numeric_or_throw<double>(value)};
    // NOLINTNEXTLINE
  } catch (...) {
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

[[nodiscard]] static toml::table normalize_attribute_map(const hic::File& hf,
                                                         const std::string& uri) {
  toml::table attributes;

  if (!uri.empty()) {
    emplace_if_valid("uri", uri, attributes);
  }

  emplace_if_valid("format", "HIC", attributes);
  emplace_if_valid("format-version", hf.version(), attributes);
  emplace_if_valid("assembly", hf.assembly(), attributes);
  emplace_if_valid("format-url", "https://github.com/aidenlab/hic-format", attributes);
  emplace_if_valid("nchroms", hf.chromosomes().remove_ALL().size(), attributes);

  for (const auto& [k, v] : hf.attributes()) {
    std::visit([&, key = k](const auto& value) { emplace_if_valid(key, value, attributes); },
               try_parse_str(v));
  }

  return attributes;
}

[[nodiscard]] static toml::table extract_top_lvl_metadata_hic(const std::filesystem::path& p,
                                                              bool include_file_path) {
  const auto resolution = hic::utils::list_resolutions(p).back();
  return normalize_attribute_map(
      hic::File(p.string(), resolution, hic::MatrixType::observed, hic::MatrixUnit::BP, 1),
      include_file_path ? p.string() : "");
}

[[nodiscard]] static toml::array read_hic_matrix_types(const std::filesystem::path& p,
                                                       std::uint32_t resolution) {
  toml::array buff;
  for (const auto& mt : {"observed", "expected", "oe"}) {
    try {
      std::ignore =
          hic::File(p.string(), resolution, hic::ParseMatrixTypeStr(mt), hic::MatrixUnit::BP, 1);
      buff.push_back(mt);
      // NOLINTNEXTLINE
    } catch (...) {
    }
  }
  return buff;
}

[[nodiscard]] static toml::array read_hic_normalizations_ev(const hic::File& hf) {
  toml::array buff;
  for (const auto& norm : hf.avail_normalizations()) {
    try {
      std::ignore = hf.expected_values(hf.chromosomes().longest_chromosome(), norm);
      buff.push_back(norm.to_string());
      // NOLINTNEXTLINE
    } catch (...) {
    }
  }

  return buff;
}

[[nodiscard]] static toml::array read_hic_normalizations(const hic::File& hf) {
  toml::array buff;
  for (const auto& norm : hf.avail_normalizations()) {
    buff.push_back(norm.to_string());
  }

  return buff;
}

[[nodiscard]] static std::vector<std::pair<std::string, toml::table>> extract_nested_metadata_hic(
    const std::filesystem::path& p) {
  std::vector<std::pair<std::string, toml::table>> nested_attributes{};
  for (const auto& resolution : hic::utils::list_resolutions(p)) {
    const hic::File hf(p.string(), resolution, hic::MatrixType::observed, hic::MatrixUnit::BP, 1);

    toml::table attributes;

    emplace_if_valid("matrix-types", read_hic_matrix_types(p, resolution), attributes);
    emplace_if_valid("nbins", hf.bins().size(), attributes);
    emplace_if_valid("normalizations", read_hic_normalizations(hf), attributes);
    emplace_if_valid("normalizations-ev", read_hic_normalizations_ev(hf), attributes);

    nested_attributes.emplace_back(fmt::to_string(resolution), std::move(attributes));
  }

  return nested_attributes;
}

int print_hic_metadata(const std::filesystem::path& p, MetadataOutputFormat format,
                       bool include_file_path, bool recursive) {
  auto attributes = extract_top_lvl_metadata_hic(p, include_file_path);
  std::vector<std::pair<std::string, toml::table>> nested_attributes{};

  toml::array resolutions;
  for (const auto& resolution : hic::utils::list_resolutions(p)) {
    resolutions.push_back(static_cast<std::int64_t>(resolution));
  }
  emplace_if_valid("resolutions", resolutions, attributes);

  if (recursive) {
    nested_attributes = extract_nested_metadata_hic(p);
  }

  print_attributes(attributes, nested_attributes, format);
  return 0;
}

}  // namespace hictk::tools
