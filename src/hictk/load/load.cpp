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

template <typename N>
void merge_coolers(std::vector<CoolerChunk<N>>& coolers, std::string_view dest, bool force) {
  if (force) {
    std::filesystem::remove(dest);
  }

  if (coolers.size() == 1) {
    spdlog::info(FMT_STRING("Moving temporary file to {}..."), dest);
    const auto path = coolers.front().clr.uri();
    coolers.front().clr.close();
    std::filesystem::copy_file(path, dest);
    return;
  }

  spdlog::info(FMT_STRING("Merging {} intermediate files into {}..."), coolers.size(), dest);
  const auto chroms = coolers.front().clr.chromosomes();
  const auto bin_size = coolers.front().clr.bin_size();

  using PixelIt = decltype(coolers.front().first);
  std::vector<PixelIt> heads;
  std::vector<PixelIt> tails;

  heads.reserve(coolers.size());
  tails.reserve(coolers.size());
  for (auto&& clr : coolers) {
    heads.emplace_back(std::move(clr.first));
    tails.emplace_back(std::move(clr.last));
  }
  coolers.clear();

  cooler::utils::merge(heads, tails, chroms, bin_size, dest, force);
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
        std::vector<CoolerChunk<N>> chunks{};
        buffer.clear();
        for (std::size_t i = 0; true; ++i) {
          const auto tmp_uri = tmpdir() / fmt::format(FMT_STRING("chunk_{:03d}.cool"), i);
          spdlog::info(FMT_STRING("writing chunk #{} to intermediate file {}..."), i + 1, tmp_uri);
          chunks.emplace_back(ingest_pixels_unsorted(
              cooler::File::create_new_cooler<N>(tmp_uri.string(), chroms, c.bin_size, c.force),
              buffer, format, c.validate_pixels));

          if (chunks.back().first == chunks.back().last) {
            chunks.pop_back();
            break;
          }
          spdlog::info(FMT_STRING("done writing to file {}..."), tmp_uri);
        }
        merge_coolers(chunks, c.uri, c.force);
      },
      write_buffer);
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
        std::vector<CoolerChunk<N>> chunks{};

        for (std::size_t i = 0; true; ++i) {
          const auto tmp_uri = tmpdir() / fmt::format(FMT_STRING("chunk_{:03d}.cool"), i);
          chunks.emplace_back(ingest_pairs_unsorted(
              cooler::File::create_new_cooler<N>(tmp_uri.string(), chroms, c.bin_size, c.force),
              buffer, c.batch_size, format, c.validate_pixels));

          if (chunks.back().first == chunks.back().last) {
            chunks.pop_back();
            break;
          }
          spdlog::info(FMT_STRING("Done writing to tmp file {}..."), tmp_uri);
        }
        merge_coolers(chunks, c.uri, c.force);
      },
      write_buffer);
}

int load_subcmd(const LoadConfig& c) {
  const auto format = format_from_string(c.format);
  const auto pixel_has_count = format == Format::COO || format == Format::BG2;

  if (c.assume_sorted && pixel_has_count) {
    spdlog::info(FMT_STRING("begin loading presorted pixels..."));
    ingest_pixels_sorted(c);
  } else if (!c.assume_sorted && pixel_has_count) {
    spdlog::info(FMT_STRING("begin loading un-sorted pixels..."));
    ingest_pixels_unsorted(c);
  } else if (c.assume_sorted && !pixel_has_count) {
    spdlog::info(FMT_STRING("begin loading presorted pairs..."));
    ingest_pairs_sorted(c);
  } else {
    spdlog::info(FMT_STRING("begin loading un-sorted pairs..."));
    ingest_pairs_unsorted(c);
  }
  return 0;
}

}  // namespace hictk::tools
