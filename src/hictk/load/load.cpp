// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>
#include <fmt/std.h>
#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/pixel.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/tools.hpp"

namespace hictk::tools {

enum class Format { COO, BG2, VP, _4DN };
[[nodiscard]] static Format format_from_string(std::string_view s) {
  if (s == "coo") {
    return Format::COO;
  }
  if (s == "bg2") {
    return Format::BG2;
  }
  if (s == "validpairs") {
    return Format::VP;
  }
  assert(s == "4dn");
  return Format::_4DN;
}

template <typename N>
[[nodiscard]] static Pixel<N> parse_pixel(const BinTable& bins, std::string_view line,
                                          Format format) {
  switch (format) {
    case Format::COO:
      return Pixel<N>::from_coo(bins, line);
    case Format::BG2:
      return Pixel<N>::from_bg2(bins, line);
    case Format::VP:
      return Pixel<N>::from_validpair(bins, line);
    case Format::_4DN:
      return Pixel<N>::from_4dn_pairs(bins, line);
  }
  HICTK_UNREACHABLE_CODE;
}

[[nodiscard]] bool line_is_header(std::string_view line) {
  return !line.empty() && line.front() == '#';
}

template <typename N>
static void read_batch(const BinTable& bins, std::vector<Pixel<N>>& buffer, Format format) {
  buffer.clear();
  std::string line{};
  try {
    while (std::getline(std::cin, line)) {
      if (line_is_header(line)) {
        continue;
      }
      if (buffer.size() == buffer.capacity()) {
        return;
      }
      buffer.emplace_back(parse_pixel<N>(bins, line, format));
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("encountered error while processing the following line:\n"
                               "\"{}\"\n"
                               "Cause: {}"),
                    line, e.what()));
  }
}

template <typename N>
static void read_and_sort_batch(const BinTable& bins, std::vector<Pixel<N>>& buffer,
                                Format format) {
  read_batch(bins, buffer, format);
  std::sort(buffer.begin(), buffer.end());
}

template <typename N>
static void read_and_aggregate_batch(const BinTable& bins, std::vector<Pixel<N>>& buffer,
                                     Format format) {
  buffer.clear();
  std::string line{};
  try {
    while (std::getline(std::cin, line)) {
      if (line_is_header(line)) {
        continue;
      }
      if (buffer.size() == buffer.capacity()) {
        return;
      }

      auto pixel = parse_pixel<N>(bins, line, format);
      if (!buffer.empty() && buffer.back().coords == pixel.coords) {
        buffer.back().count += pixel.count;
      } else {
        buffer.emplace_back(std::move(pixel));
      }
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("encountered error while processing the following line:\n"
                               "\"{}\"\n"
                               "Cause: {}"),
                    line, e.what()));
  }
}

template <typename N>
static void read_sort_and_aggregate_batch(const BinTable& bins, std::vector<Pixel<N>>& flat_buffer,
                                          Format format) {
  phmap::btree_map<PixelCoordinates, N> buffer;

  std::string line{};
  try {
    while (std::getline(std::cin, line)) {
      if (line_is_header(line)) {
        continue;
      }
      if (buffer.size() == flat_buffer.capacity()) {
        return;
      }
      auto pixel = parse_pixel<N>(bins, line, format);
      auto [node, inserted] = buffer.try_emplace(std::move(pixel.coords), pixel.count);
      if (!inserted) {
        node->second += pixel.count;
      }
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("encountered error while processing the following line:\n"
                               "\"{}\"\n"
                               "Cause: {}"),
                    line, e.what()));
  }

  flat_buffer.resize(buffer.size());
  std::transform(std::make_move_iterator(buffer.begin()), std::make_move_iterator(buffer.end()),
                 flat_buffer.begin(), [](auto&& pair) {
                   return Pixel<N>{std::move(pair.first), pair.second};
                 });
}

template <typename N>
static void ingest_pixels_sorted(hictk::cooler::File&& clr, Format format, std::size_t batch_size) {
  std::vector<Pixel<N>> buffer(batch_size);

  std::string line{};
  std::size_t i = 0;
  try {
    while (!std::cin.eof()) {
      read_batch(clr.bins(), buffer, format);
      clr.append_pixels(buffer.begin(), buffer.end(), true);
      buffer.clear();
    }
    if (!buffer.empty()) {
      clr.append_pixels(buffer.begin(), buffer.end(), true);
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while processing line(s) {}-{}: {}"),
                    i - buffer.size(), i, e.what()));
  }
}

template <typename N>
static std::string ingest_pixels_unsorted(hictk::cooler::File&& clr,
                                          std::vector<Pixel<N>>& write_buffer, Format format) {
  assert(write_buffer.capacity() != 0);

  read_and_sort_batch(clr.bins(), write_buffer, format);

  if (write_buffer.empty()) {
    assert(std::cin.eof());
    return "";
  }

  clr.append_pixels(write_buffer.begin(), write_buffer.end(), true);
  write_buffer.clear();

  return clr.uri();
}

template <typename N>
static void ingest_pairs_sorted(hictk::cooler::File&& clr, Format format, std::size_t batch_size) {
  std::vector<Pixel<N>> buffer(batch_size);

  std::size_t i = 0;
  try {
    while (!std::cin.eof()) {
      read_and_aggregate_batch(clr.bins(), buffer, format);
      clr.append_pixels(buffer.begin(), buffer.end(), true);
      buffer.clear();
    }
    if (!buffer.empty()) {
      clr.append_pixels(buffer.begin(), buffer.end(), true);
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while processing line(s) {}-{}: {}"),
                    i - buffer.size(), i, e.what()));
  }
}

template <typename N>
static std::string ingest_pairs_unsorted(hictk::cooler::File&& clr,
                                         std::vector<Pixel<N>>& write_buffer, Format format) {
  assert(write_buffer.capacity() != 0);

  read_sort_and_aggregate_batch(clr.bins(), write_buffer, format);

  if (write_buffer.empty()) {
    assert(std::cin.eof());
    return "";
  }

  clr.append_pixels(write_buffer.begin(), write_buffer.end());
  write_buffer.clear();

  return clr.uri();
}

void merge_coolers(const std::vector<std::string>& sources, std::string_view dest, bool force,
                   std::uint32_t verbosity) {
  if (force) {
    std::filesystem::remove(dest);
  }

  if (sources.size() == 1) {
    spdlog::info(FMT_STRING("Moving temporary file to {}..."), dest);
    std::filesystem::copy_file(sources.front(), dest);
    return;
  }

  spdlog::info(FMT_STRING("Merging {} intermediate files into {}..."), sources.size(), dest);
  hictk::cooler::utils::merge(sources.begin(), sources.end(), dest, force, 500'000, verbosity != 3);
}

void ingest_pixels_sorted(const LoadConfig& c) {
  assert(c.assume_sorted);
  auto chroms = Reference(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  c.count_as_float
      ? ingest_pixels_sorted<double>(
            hictk::cooler::File::create_new_cooler<double>(c.uri, chroms, c.bin_size, c.force),
            format, c.batch_size)
      : ingest_pixels_sorted<std::int32_t>(hictk::cooler::File::create_new_cooler<std::int32_t>(
                                               c.uri, chroms, c.bin_size, c.force),
                                           format, c.batch_size);
}

void ingest_pixels_unsorted(const LoadConfig& c) {
  assert(!c.assume_sorted);
  auto chroms = Reference(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  const internal::TmpDir tmpdir{};
  std::vector<std::string> uris{};

  using IntBuff = std::vector<Pixel<std::int32_t>>;
  using FPBuff = std::vector<Pixel<double>>;
  std::variant<IntBuff, FPBuff> write_buffer{};
  if (c.count_as_float) {
    write_buffer = FPBuff(c.batch_size);
  } else {
    write_buffer = IntBuff(c.batch_size);
  }

  std::visit(
      [&](auto& buffer) {
        using N = decltype(buffer.front().count);
        for (std::size_t i = 0; true; ++i) {
          const auto tmp_uri = tmpdir() / fmt::format(FMT_STRING("chunk_{:03d}.cool"), i);
          uris.emplace_back(
              ingest_pixels_unsorted(hictk::cooler::File::create_new_cooler<N>(
                                         tmp_uri.string(), chroms, c.bin_size, c.force),
                                     buffer, format));
          if (uris.back().empty()) {
            uris.pop_back();
            break;
          }
          spdlog::info(FMT_STRING("Done writing to tmp file {}..."), tmp_uri);
        }
      },
      write_buffer);

  merge_coolers(uris, c.uri, c.force, c.verbosity);
}

void ingest_pairs_sorted(const LoadConfig& c) {
  assert(c.assume_sorted);
  auto chroms = Reference(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  c.count_as_float
      ? ingest_pairs_sorted<double>(
            hictk::cooler::File::create_new_cooler<double>(c.uri, chroms, c.bin_size, c.force),
            format, c.batch_size)
      : ingest_pairs_sorted<std::int32_t>(hictk::cooler::File::create_new_cooler<std::int32_t>(
                                              c.uri, chroms, c.bin_size, c.force),
                                          format, c.batch_size);
}

void ingest_pairs_unsorted(const LoadConfig& c) {
  assert(!c.assume_sorted);
  auto chroms = Reference(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  const internal::TmpDir tmpdir{};
  std::vector<std::string> uris{};

  using IntBuff = std::vector<Pixel<std::int32_t>>;
  using FPBuff = std::vector<Pixel<double>>;
  std::variant<IntBuff, FPBuff> write_buffer{};
  if (c.count_as_float) {
    write_buffer = FPBuff(c.batch_size);
  } else {
    write_buffer = IntBuff(c.batch_size);
  }

  std::visit(
      [&](auto& buffer) {
        using N = decltype(buffer.front().count);
        for (std::size_t i = 0; true; ++i) {
          const auto tmp_uri = tmpdir() / fmt::format(FMT_STRING("chunk_{:03d}.cool"), i);
          uris.emplace_back(
              ingest_pairs_unsorted(hictk::cooler::File::create_new_cooler<N>(
                                        tmp_uri.string(), chroms, c.bin_size, c.force),
                                    buffer, format));
          if (uris.back().empty()) {
            uris.pop_back();
            break;
          }
          spdlog::info(FMT_STRING("Done writing to tmp file {}..."), tmp_uri);
        }
      },
      write_buffer);

  merge_coolers(uris, c.uri, c.force, c.verbosity);
}

void load_subcmd(const LoadConfig& c) {
  const auto format = format_from_string(c.format);
  const auto pixel_has_count = format == Format::COO || format == Format::BG2;

  if (c.assume_sorted && pixel_has_count) {
    ingest_pixels_sorted(c);
  } else if (!c.assume_sorted && pixel_has_count) {
    ingest_pixels_unsorted(c);
  } else if (c.assume_sorted && !pixel_has_count) {
    ingest_pairs_sorted(c);
  } else {
    ingest_pairs_unsorted(c);
  }
}

}  // namespace hictk::tools
