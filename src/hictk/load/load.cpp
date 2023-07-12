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

#include "./common.hpp"
#include "./load_pairs.hpp"
#include "./load_pixels.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/pixel.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/tools.hpp"

namespace hictk::tools {

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
  cooler::utils::merge(sources.begin(), sources.end(), dest, force, 500'000, verbosity != 3);
}

void ingest_pixels_sorted(const LoadConfig& c) {
  assert(c.assume_sorted);
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  c.count_as_float
      ? ingest_pixels_sorted<double>(
            cooler::File::create_new_cooler<double>(c.uri, chroms, c.bin_size, c.force), format,
            c.batch_size, c.validate_pixels)
      : ingest_pixels_sorted<std::int32_t>(
            cooler::File::create_new_cooler<std::int32_t>(c.uri, chroms, c.bin_size, c.force),
            format, c.batch_size, c.validate_pixels);
}

void ingest_pixels_unsorted(const LoadConfig& c) {
  assert(!c.assume_sorted);
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  const internal::TmpDir tmpdir{};
  std::vector<std::string> uris{};

  using IntBuff = std::vector<ThinPixel<std::int32_t>>;
  using FPBuff = std::vector<ThinPixel<double>>;
  std::variant<IntBuff, FPBuff> write_buffer{};
  if (c.count_as_float) {
    write_buffer = FPBuff(c.batch_size);
  } else {
    write_buffer = IntBuff(c.batch_size);
  }

  std::visit(
      [&](auto& buffer) {
        using N = decltype(buffer.front().count);
        buffer.clear();
        for (std::size_t i = 0; true; ++i) {
          const auto tmp_uri = tmpdir() / fmt::format(FMT_STRING("chunk_{:03d}.cool"), i);
          uris.emplace_back(ingest_pixels_unsorted(
              cooler::File::create_new_cooler<N>(tmp_uri.string(), chroms, c.bin_size, c.force),
              buffer, format, c.validate_pixels));
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
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  c.count_as_float
      ? ingest_pairs_sorted<double>(
            cooler::File::create_new_cooler<double>(c.uri, chroms, c.bin_size, c.force), format,
            c.batch_size, c.validate_pixels)
      : ingest_pairs_sorted<std::int32_t>(
            cooler::File::create_new_cooler<std::int32_t>(c.uri, chroms, c.bin_size, c.force),
            format, c.batch_size, c.validate_pixels);
}

static void ingest_pairs_unsorted(const LoadConfig& c) {
  assert(!c.assume_sorted);
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  const internal::TmpDir tmpdir{};
  std::vector<std::string> uris{};

  using IntBuff = std::vector<ThinPixel<std::int32_t>>;
  using FPBuff = std::vector<ThinPixel<double>>;
  std::variant<IntBuff, FPBuff> write_buffer{};
  if (c.count_as_float) {
    write_buffer = FPBuff{};
  } else {
    write_buffer = IntBuff{};
  }

  std::visit(
      [&](auto& buffer) {
        using N = decltype(buffer.begin()->count);

        for (std::size_t i = 0; true; ++i) {
          const auto tmp_uri = tmpdir() / fmt::format(FMT_STRING("chunk_{:03d}.cool"), i);
          uris.emplace_back(ingest_pairs_unsorted(
              cooler::File::create_new_cooler<N>(tmp_uri.string(), chroms, c.bin_size, c.force),
              buffer, c.batch_size, format, c.validate_pixels));
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

int load_subcmd(const LoadConfig& c) {
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
  return 1;
}

}  // namespace hictk::tools
