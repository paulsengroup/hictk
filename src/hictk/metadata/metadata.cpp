// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/string_utils.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/type_traits.hpp"
#include "hictk/version.hpp"

namespace hictk::tools {

enum class MetadataOutputFormat : std::uint8_t { json, toml, tsv, yaml };

class Attribute {
  std::string _key{};
  std::optional<std::string> _value{};
  bool _quote{};

 public:
  Attribute() = default;
  template <typename T>
  Attribute(std::string_view key, const T& value) : _key(key) {
    auto [str, needs_quoting] = to_str(value);
    _value = std::move(str);  // NOLINT
    _quote = needs_quoting;   // NOLINT
  }

  explicit operator bool() const noexcept { return !_key.empty() && _value.has_value(); }

  bool operator<(const Attribute& other) const noexcept {
    if (_key == other._key) {
      return _value < other._value;
    }
    return _key < other._key;
  }

  template <typename T>
  [[nodiscard]] static std::pair<std::optional<std::string>, bool> to_str(const T& value) {
    return std::make_pair(std::make_optional(fmt::to_string(value)), false);
  }

  [[nodiscard]] static std::pair<std::optional<std::string>, bool> to_str(
      const std::string& value) {
    try {
      // NOLINTNEXTLINE
      return to_str(internal::parse_numeric_or_throw<std::uint64_t>(value));
    } catch (...) {  // NOLINT
    }

    try {
      // NOLINTNEXTLINE
      return to_str(internal::parse_numeric_or_throw<std::int64_t>(value));
    } catch (...) {  // NOLINT
    }

    try {
      // NOLINTNEXTLINE
      return to_str(internal::parse_numeric_or_throw<double>(value));
    } catch (...) {  // NOLINT
    }

    if (value == "true" || value == "True") {
      return to_str(true);
    }
    if (value == "false" || value == "False") {
      return to_str(false);
    }

    if (value == "NULL" || value == "Null" || value == "null" || value == "None") {
      return std::make_pair(std::make_optional("null"), false);
    }

    return std::make_pair(std::make_optional(value), true);
  }

  template <typename... T>
  [[nodiscard]] static std::pair<std::optional<std::string>, bool> to_str(
      const std::variant<T...>& var) {
    return std::visit([&](const auto& value) { return to_str(value); }, var);
  }

  template <typename T>
  [[nodiscard]] static std::pair<std::optional<std::string>, bool> to_str(
      const std::optional<T>& opt) {
    if (opt.has_value()) {
      return to_str(opt.value());
    }
    return {};
  }

  template <MetadataOutputFormat out_fmt,
            typename std::enable_if_t<out_fmt == MetadataOutputFormat::json>* = nullptr>
  [[nodiscard]] std::string format() const {
    if (!_value.has_value()) {
      return "";
    }

    if (*_value == "{}") {
      return fmt::format(FMT_STRING("\t{:?}: {{}}"), _key);
    }

    std::string str{};
    if (_value->find('\n') != std::string_view::npos) {
      // indent value
      str = internal::str_replace(*_value, "\n", "\t\n");
      str = std::string{internal::remove_suffix(str, "\t\n")};
    } else {
      str = *_value;
    }

    return _quote ? fmt::format(FMT_STRING("\t{:?}: \"{}\""), _key, internal::escape_str(str))
                  : fmt::format(FMT_STRING("\t{:?}: {}"), _key, internal::escape_str(str));
  }

  template <MetadataOutputFormat out_fmt,
            typename std::enable_if_t<out_fmt == MetadataOutputFormat::toml>* = nullptr>
  [[nodiscard]] std::string format() const {
    if (!_value.has_value() || (_key == "metadata" && *_value == "{}")) {
      return "";
    }

    auto str = internal::escape_str(*_value);

    if (_quote && str.find("\\n") != std::string::npos) {
      // Format multi-line str
      return fmt::format(FMT_STRING("{} = \"\"\"\n{}\"\"\""), _key, str);
    }

    return _quote ? fmt::format(FMT_STRING("{} = \"{}\""), _key, str)
                  : fmt::format(FMT_STRING("{} = {}"), _key, str);
  }

  template <MetadataOutputFormat out_fmt,
            typename std::enable_if_t<out_fmt == MetadataOutputFormat::tsv>* = nullptr>
  [[nodiscard]] std::string format() const {
    if (!_value.has_value()) {
      return "";
    }

    if (_key == "metadata" && *_value == "{}") {
      return fmt::format(FMT_STRING("{:?}\t"), _key);
    }

    if (_quote) {
      return fmt::format(FMT_STRING("{:?}\t{:?}"), _key, *_value);
    }
    return fmt::format(FMT_STRING("{:?}\t{}"), _key, *_value);
  }

  template <MetadataOutputFormat out_fmt,
            typename std::enable_if_t<out_fmt == MetadataOutputFormat::yaml>* = nullptr>
  [[nodiscard]] std::string format() const {
    if (!_value.has_value()) {
      return "";
    }

    if (_key == "metadata" && *_value == "{}") {
      return fmt::format(FMT_STRING("{}: null"), _key);
    }

    auto str = internal::escape_str(*_value);
    str = internal::str_replace(str, "\\r\\n", "\t\n");
    str = internal::str_replace(str, "\\n", "\t\n");

    if (str.find('\n') == std::string::npos) {
      return fmt::format(FMT_STRING("{}: {}"), _key, str);
    }
    return fmt::format(FMT_STRING("{}: |\n  {}"), _key, str);
  }
};

[[nodiscard]] static MetadataOutputFormat parse_output_format(std::string_view format) {
  if (format == "json") {
    return MetadataOutputFormat::json;
  }
  if (format == "toml") {
    return MetadataOutputFormat::toml;
  }
  if (format == "tsv") {
    return MetadataOutputFormat::tsv;
  }
  assert(format == "yaml");
  return MetadataOutputFormat::yaml;
}

[[nodiscard]] static std::vector<Attribute> normalize_attribute_map(
    const phmap::flat_hash_map<std::string, std::string>& map) {
  std::vector<Attribute> attributes(map.size());
  attributes.clear();

  for (const auto& [k, v] : map) {
    attributes.emplace_back(k, v);
  }

  attributes.erase(
      std::remove_if(attributes.begin(), attributes.end(), [](const auto& attr) { return !attr; }),
      attributes.end());
  std::sort(attributes.begin(), attributes.end());

  return attributes;
}

static void emplace_if_valid(Attribute&& attribute, std::vector<Attribute>& buff) {
  if (!!attribute) {
    buff.emplace_back(std::move(attribute));
  }
}

[[nodiscard]] static std::vector<Attribute> normalize_attribute_map(
    const cooler::MultiResAttributes& map) {
  std::vector<Attribute> attributes(3);
  attributes.clear();

  emplace_if_valid(Attribute("bin-type", map.bin_type), attributes);
  emplace_if_valid(Attribute("format", map.format), attributes);
  emplace_if_valid(Attribute("format-version", map.format_version), attributes);

  std::sort(attributes.begin(), attributes.end());

  return attributes;
}

[[nodiscard]] static std::vector<Attribute> normalize_attribute_map(
    const cooler::SingleCellAttributes& map) {
  std::vector<Attribute> attributes(13);
  attributes.clear();

  emplace_if_valid(Attribute("bin-size", map.bin_size), attributes);
  emplace_if_valid(Attribute("bin-type", map.bin_type), attributes);
  emplace_if_valid(Attribute("format", map.format), attributes);
  emplace_if_valid(Attribute("format-version", map.format_version), attributes);

  emplace_if_valid(Attribute("creation-date", map.creation_date), attributes);
  emplace_if_valid(Attribute("generated-by", map.generated_by), attributes);
  emplace_if_valid(Attribute("assembly", map.assembly), attributes);
  emplace_if_valid(Attribute("metadata", map.metadata), attributes);

  emplace_if_valid(Attribute("format-url", map.format_url), attributes);
  emplace_if_valid(Attribute("nbins", map.nbins), attributes);
  emplace_if_valid(Attribute("ncells", map.ncells), attributes);
  emplace_if_valid(Attribute("nchroms", map.nchroms), attributes);
  emplace_if_valid(Attribute("storage-mode", map.storage_mode), attributes);

  std::sort(attributes.begin(), attributes.end());

  return attributes;
}

[[nodiscard]] static std::vector<Attribute> normalize_attribute_map(const cooler::Attributes& map) {
  std::vector<Attribute> attributes(15);
  attributes.clear();

  emplace_if_valid(Attribute("bin-size", map.bin_size), attributes);
  emplace_if_valid(Attribute("bin-type", map.bin_type), attributes);
  emplace_if_valid(Attribute("format", map.format), attributes);
  emplace_if_valid(Attribute("format-version", map.format_version), attributes);
  emplace_if_valid(Attribute("storage-mode", map.storage_mode), attributes);

  emplace_if_valid(Attribute("creation-date", map.creation_date), attributes);
  emplace_if_valid(Attribute("generated-by", map.generated_by), attributes);
  emplace_if_valid(Attribute("assembly", map.assembly), attributes);
  emplace_if_valid(Attribute("metadata", map.metadata), attributes);

  emplace_if_valid(Attribute("format-url", map.format_url), attributes);
  emplace_if_valid(Attribute("nbins", map.nbins), attributes);
  emplace_if_valid(Attribute("nchroms", map.nchroms), attributes);
  emplace_if_valid(Attribute("nnz", map.nnz), attributes);
  emplace_if_valid(Attribute("sum", map.sum), attributes);
  emplace_if_valid(Attribute("cis", map.cis), attributes);

  std::sort(attributes.begin(), attributes.end());

  return attributes;
}

static void format_to_json(const std::vector<Attribute>& attrs) {
  fmt::print(FMT_STRING("{{\n"));

  for (std::size_t i = 0; i < attrs.size() - 1; ++i) {
    if (const auto s = attrs[i].format<MetadataOutputFormat::json>(); !s.empty()) {
      fmt::print(FMT_STRING("{},\n"), s);
    }
  }
  if (const auto s = attrs.back().format<MetadataOutputFormat::json>(); !s.empty()) {
    fmt::print(FMT_STRING("{}\n}}"), s);
  }
}

static void format_to_tsv(const std::vector<Attribute>& attrs) {
  fmt::print(FMT_STRING("\"attribute\"\t\"value\"\n"));

  for (const auto& a : attrs) {
    if (const auto s = a.format<MetadataOutputFormat::tsv>(); !s.empty()) {
      fmt::print(FMT_STRING("{}\n"), s);
    }
  }
}

static void format_to_toml(const std::vector<Attribute>& attrs) {
  fmt::print(FMT_STRING("# Metadata generated by {}\n"), hictk::config::version::str_long());

  for (const auto& a : attrs) {
    if (const auto s = a.format<MetadataOutputFormat::toml>(); !s.empty()) {
      fmt::print(FMT_STRING("{}\n"), s);
    }
  }
}

static void format_to_yaml(const std::vector<Attribute>& attrs) {
  fmt::print(FMT_STRING("--- # Metadata generated by {}\n"), hictk::config::version::str_long());

  for (const auto& a : attrs) {
    if (const auto s = a.format<MetadataOutputFormat::yaml>(); !s.empty()) {
      fmt::print(FMT_STRING("{}\n"), s);
    }
  }
}

static void print_attributes(const std::vector<Attribute>& attrs, MetadataOutputFormat format) {
  switch (format) {
    case MetadataOutputFormat::json:
      return format_to_json(attrs);  // NOLINT
    case MetadataOutputFormat::toml:
      return format_to_toml(attrs);  // NOLINT
    case MetadataOutputFormat::tsv:
      return format_to_tsv(attrs);  // NOLINT
    case MetadataOutputFormat::yaml:
      return format_to_yaml(attrs);  // NOLINT
  }
}

[[nodiscard]] static int print_hic_metadata(const std::filesystem::path& p,
                                            MetadataOutputFormat format) {
  const auto resolution = hic::utils::list_resolutions(p).front();
  const auto attributes = normalize_attribute_map(hic::File(p, resolution).attributes());
  print_attributes(attributes, format);

  return 0;
}

[[nodiscard]] static int print_mcool_metadata(const std::filesystem::path& p,
                                              MetadataOutputFormat format) {
  const auto attributes = normalize_attribute_map(cooler::MultiResFile(p).attributes());
  print_attributes(attributes, format);

  return 0;
}

[[nodiscard]] static int print_scool_metadata(const std::filesystem::path& p,
                                              MetadataOutputFormat format) {
  const auto attributes = normalize_attribute_map(cooler::SingleCellFile(p).attributes());
  print_attributes(attributes, format);

  return 0;
}

[[nodiscard]] static int print_cool_metadata(const std::filesystem::path& p,
                                             MetadataOutputFormat format) {
  const auto attributes = normalize_attribute_map(cooler::File(p.string()).attributes());
  print_attributes(attributes, format);

  return 0;
}

int metadata_subcmd(const MetadataConfig& c) {
  const auto output_format = parse_output_format(c.output_format);
  if (c.input_format == "hic") {
    return print_hic_metadata(c.uri, output_format);
  }
  if (c.input_format == "mcool") {
    return print_mcool_metadata(c.uri, output_format);
  }
  if (c.input_format == "scool") {
    return print_scool_metadata(c.uri, output_format);
  }
  assert(c.input_format == "cool");
  return print_cool_metadata(c.uri, output_format);
}

}  // namespace hictk::tools
