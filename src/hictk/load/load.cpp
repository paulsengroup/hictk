// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <variant>
#include <vector>

#include "./common.hpp"
#include "./load_cooler.hpp"
#include "./load_hic.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/tools.hpp"

namespace hictk::tools {

[[nodiscard]] static BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes,
                                             std::uint32_t bin_size) {
  auto chroms = Reference::from_chrom_sizes(path_to_chrom_sizes);
  return {chroms, bin_size};
}

[[nodiscard]] static BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes,
                                             const std::filesystem::path& path_to_bin_table) {
  auto chroms = Reference::from_chrom_sizes(path_to_chrom_sizes);

  std::ifstream ifs{};
  ifs.exceptions(std::ios::badbit);
  ifs.open(path_to_bin_table);

  std::vector<std::uint32_t> start_pos{};
  std::vector<std::uint32_t> end_pos{};

  std::string line{};
  GenomicInterval record{};
  bool fixed_bin_size = true;
  std::uint32_t bin_size = 0;
  while (std::getline(ifs, line)) {
    record = GenomicInterval::parse_bed(chroms, line);
    if (bin_size == 0) {
      bin_size = record.size();
    }

    fixed_bin_size &= record.size() == bin_size || record.chrom().size() == record.end();

    start_pos.push_back(record.start());
    end_pos.push_back(record.end());
  }

  if (fixed_bin_size) {
    SPDLOG_INFO(FMT_STRING("detected bin table with uniform bin size."));
    return {chroms, bin_size};
  }

  SPDLOG_INFO(FMT_STRING("detected bin table with variable bin size."));
  return {chroms, start_pos, end_pos};
}

static Stats ingest_pixels_hic(const LoadConfig& c) {
  const auto format = format_from_string(c.format);
  const auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  return ingest_pixels_hic(c.output_path, c.tmp_dir, chroms, c.bin_size, c.assembly, c.offset,
                           format, c.threads, c.batch_size, c.compression_lvl, c.force);
}

static Stats ingest_pixels_cooler(const LoadConfig& c) {
  assert(c.output_format == "cool");
  const auto format = format_from_string(c.format);
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto tmp_cooler_path =
      (c.tmp_dir / (std::filesystem::path{c.output_path}.filename().string() + ".tmp")).string();

  return c.assume_sorted ? ingest_pixels_sorted_cooler(c.output_path, chroms, c.bin_size, c.offset,
                                                       format, c.batch_size, c.force,
                                                       c.count_as_float, c.validate_pixels)
                         : ingest_pixels_unsorted_cooler(
                               c.output_path, tmp_cooler_path, chroms, c.bin_size, c.offset, format,
                               c.batch_size, c.force, c.count_as_float, c.validate_pixels);
}

static Stats ingest_pairs_cooler(const LoadConfig& c) {
  auto bins = c.path_to_bin_table.empty()
                  ? init_bin_table(c.path_to_chrom_sizes, c.bin_size)
                  : init_bin_table(c.path_to_chrom_sizes, c.path_to_bin_table);
  const auto format = format_from_string(c.format);
  const auto tmp_cooler_path =
      (c.tmp_dir / (std::filesystem::path{c.output_path}.filename().string() + ".tmp")).string();

  return ingest_pairs_cooler(c.output_path, tmp_cooler_path, bins, c.offset, format, c.batch_size,
                             c.force, c.count_as_float, c.validate_pixels);
}

static Stats ingest_pairs_hic(const LoadConfig& c) {
  const auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  return ingest_pairs_hic(c.output_path, c.tmp_dir, chroms, c.bin_size, c.assembly, c.offset,
                          format, c.threads, c.batch_size, c.compression_lvl, c.force);
}

static Stats ingest_pixels(const LoadConfig& c) {
  if (c.output_format == "hic") {
    return ingest_pixels_hic(c);
  }

  return ingest_pixels_cooler(c);
}

static Stats ingest_pairs(const LoadConfig& c) {
  if (c.output_format == "hic") {
    return ingest_pairs_hic(c);
  }

  return ingest_pairs_cooler(c);
}

int load_subcmd(const LoadConfig& c) {
  const auto format = format_from_string(c.format);
  const auto pixel_has_count = format == Format::COO || format == Format::BG2;
  const auto t0 = std::chrono::system_clock::now();

  const auto stats = pixel_has_count ? ingest_pixels(c) : ingest_pairs(c);

  const auto t1 = std::chrono::system_clock::now();
  const auto delta = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();

  std::visit(
      [&](const auto& sum) {
        SPDLOG_INFO(FMT_STRING("ingested {} interactions ({} nnz) in {}s!"), sum, stats.nnz,
                    static_cast<double>(delta) / 1.0e9);
      },
      stats.sum);

  return 0;
}

}  // namespace hictk::tools
