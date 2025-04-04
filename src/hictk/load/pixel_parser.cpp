// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./pixel_parser.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "./common.hpp"
#include "./init_bin_table.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/reference.hpp"
#include "hictk/string_utils.hpp"
#include "hictk/tools/compressed_io.hpp"

namespace hictk::tools {

PixelParser::PixelParser(const std::filesystem::path& path_, std::uint32_t resolution,
                         Format format_, std::string_view assembly, bool drop_unknown_chroms)
    : _reader(path_ != "-" ? io::CompressedReader{path_} : io::CompressedReader{}),
      _format(format_),
      _drop_unknown_chroms(drop_unknown_chroms) {
  auto header = parse_header(_format);
  if (!header.chromosomes) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to read chromosomes from {}"), path()));
  }

  if (header.assembly.empty() || assembly != "unknown") {
    _assembly = std::string{assembly};
  } else {
    _assembly = std::move(header.assembly);
  }

  assert(resolution != 0);
  _bins = BinTable(std::move(*header.chromosomes), resolution);
}

PixelParser::PixelParser(const std::filesystem::path& path_, BinTable bins_, Format format_,
                         std::string_view assembly_, bool drop_unknown_chroms)
    : _reader(path_ != "-" ? io::CompressedReader{path_} : io::CompressedReader{}),
      _format(format_),
      _bins(std::move(bins_)),
      _assembly(std::string{assembly_}),
      _drop_unknown_chroms(drop_unknown_chroms) {
  try {
    std::ignore = parse_header(_format);
  } catch (const std::exception& e) {
    SPDLOG_WARN(FMT_STRING("encountered an error while parsing the file header: {}"), e.what());
  }
}

const std::filesystem::path& PixelParser::path() const noexcept {
  static const std::filesystem::path stdin_{"stdin"};
  if (!!_reader) {
    return stdin_;
  }
  return _reader.path();
}

std::string_view PixelParser::assembly() const noexcept { return _assembly; }

const BinTable& PixelParser::bins() const noexcept { return _bins; }

bool PixelParser::getline(char sep) {
  const auto ok =
      !_reader ? !!std::getline(std::cin, _strbuff, sep) : _reader.getline(_strbuff, sep);
  if (!ok) {
    return ok;
  }

  if (_strbuff.empty()) {
    return getline(sep);
  }

  if (_strbuff.back() == '\r') {
    _strbuff.pop_back();
  }

  return ok;
}

std::string_view PixelParser::remove_prefix(std::string_view s, std::string_view prefix,
                                            bool strip_whitespaces) noexcept {
  if (internal::starts_with(s, prefix)) {
    s.remove_prefix(prefix.size());
  }

  while (strip_whitespaces && !s.empty()) {
    // NOLINTNEXTLINE(readability-implicit-bool-conversion)
    if (std::isspace(s.front())) {
      s.remove_prefix(1);
    } else {
      break;
    }
  }

  return s;  // NOLINT
}

std::pair<std::string, std::uint32_t> PixelParser::parse_chromsize(std::string_view line) {
  assert(internal::starts_with(line, "#chromsize:"));
  line = remove_prefix(line, "#chromsize:");
  // NOLINTNEXTLINE(readability-qualified-auto)
  auto it = std::find_if(line.begin(), line.end(), [](const char c) { return std::isspace(c); });
  if (it == line.begin() || it == line.end()) {
    throw std::runtime_error(fmt::format(FMT_STRING("malformed chromsize entry \"{}\"."), line));
  }
  std::string chrom_name{line.begin(), it};

  // NOLINTNEXTLINE(readability-implicit-bool-conversion)
  it = std::find_if(it, line.end(), [](const char c) { return !std::isspace(c); });
  const auto strlen = static_cast<std::size_t>(std::distance(it, line.end()));
  const auto offset = static_cast<std::size_t>(std::distance(line.begin(), it));

  const std::string_view chrom_size{line.data() + offset, strlen};
  if (chrom_size.empty()) {
    throw std::runtime_error(fmt::format(FMT_STRING("malformed chromsize entry \"{}\"."), line));
  }

  try {
    return std::make_pair(std::move(chrom_name),
                          internal::parse_numeric_or_throw<std::uint32_t>(chrom_size));
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("malformed chromsize entry \"{}\": {}"), line, e.what()));
  }
}

Header PixelParser::parse_header(Format format_) {
  if (format_ != Format::_4DN) {
    std::ignore = getline();
    return {};
  }

  std::optional<std::string> assembly{};
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};
  for (std::size_t i = 0; getline(); ++i) {
    try {
      if (i == 0 && !internal::starts_with(_strbuff, "## pairs format v1.0")) {
        throw std::runtime_error(
            "invalid header: first line in input file does not start with \"## pairs format "
            "v1.0\"");
      }
      if (_strbuff.empty()) {
        continue;
      }

      if (_strbuff.front() != '#') {
        break;
      }

      if (internal::starts_with(_strbuff, "#genome_assembly:")) {
        if (assembly.has_value()) {
          throw std::runtime_error("found duplicate entry for \"genome_assembly\" in file header.");
        }
        assembly = std::string{remove_prefix(_strbuff, "#genome_assembly:")};
        continue;
      }

      if (!internal::starts_with(_strbuff, "#chromsize:")) {
        continue;
      }
      auto [chrom_name, chrom_size] = parse_chromsize(_strbuff);

      chrom_names.emplace_back(std::move(chrom_name));
      chrom_sizes.push_back(chrom_size);

    } catch (const std::exception& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to parse line {} from {}: {}"), i + 1, path(), e.what()));
    }
  }

  if (!assembly.has_value()) {
    assembly = "unknown";
  }

  return {*assembly, Reference{chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()}};
}

PixelParser init_pixel_parser(Format format, const std::filesystem::path& path_to_interactions,
                              const std::filesystem::path& path_to_chrom_sizes,
                              const std::filesystem::path& path_to_bins, std::uint32_t resolution,
                              std::string_view assembly) {
  assert(format == Format::_4DN || !path_to_chrom_sizes.empty() || !path_to_bins.empty());
  BinTable bins{};

  if (!path_to_bins.empty()) {
    bins = init_bin_table(path_to_chrom_sizes, path_to_bins, resolution);
  } else if (!path_to_chrom_sizes.empty()) {
    bins = BinTable{Reference::from_chrom_sizes(path_to_chrom_sizes), resolution};
  }

  switch (format) {
    case Format::_4DN: {
      if (bins.empty()) {
        assert(resolution != 0);
        return {path_to_interactions, resolution, format, assembly};
      }
      return {path_to_interactions, std::move(bins), format, assembly};
    }
    default:
      return {path_to_interactions, std::move(bins), format, assembly};
  }
  unreachable_code();
}

}  // namespace hictk::tools
