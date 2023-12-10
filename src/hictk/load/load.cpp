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
#include "./load_pairs.hpp"
#include "./load_pixels.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/tools.hpp"

namespace hictk::tools {

void ingest_pixels_sorted(const LoadConfig& c) {
  assert(c.assume_sorted);
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  c.count_as_float ? ingest_pixels_sorted<double>(
                         cooler::File::create<double>(c.uri, chroms, c.bin_size, c.force), format,
                         c.batch_size, c.validate_pixels)
                   : ingest_pixels_sorted<std::int32_t>(
                         cooler::File::create<std::int32_t>(c.uri, chroms, c.bin_size, c.force),
                         format, c.batch_size, c.validate_pixels);
}

void ingest_pixels_unsorted(const LoadConfig& c) {
  assert(!c.assume_sorted);
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  const auto tmp_cooler_path = c.uri + ".tmp";

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
        {
          auto tmp_clr =
              cooler::SingleCellFile::create(tmp_cooler_path, chroms, c.bin_size, c.force);
          for (std::size_t i = 0; true; ++i) {
            SPDLOG_INFO(FMT_STRING("writing chunk #{} to intermediate file \"{}\"..."), i + 1,
                        tmp_cooler_path);
            const auto nnz = ingest_pixels_unsorted(tmp_clr.create_cell<N>(fmt::to_string(i)),
                                                    buffer, format, c.validate_pixels);
            SPDLOG_INFO(FMT_STRING("done writing chunk #{} to tmp file \"{}\"."), i + 1,
                        tmp_cooler_path);
            if (nnz == 0) {
              break;
            }
          }
        }
        const cooler::SingleCellFile tmp_clr(tmp_cooler_path);
        SPDLOG_INFO(FMT_STRING("merging {} chunks into \"{}\"..."), tmp_clr.cells().size(), c.uri);
        tmp_clr.aggregate<N>(c.uri, c.force);
      },
      write_buffer);
  std::filesystem::remove(tmp_cooler_path);
}

void ingest_pairs_sorted(const LoadConfig& c) {
  assert(c.assume_sorted);
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  c.count_as_float ? ingest_pairs_sorted<double>(
                         cooler::File::create<double>(c.uri, chroms, c.bin_size, c.force), format,
                         c.batch_size, c.validate_pixels)
                   : ingest_pairs_sorted<std::int32_t>(
                         cooler::File::create<std::int32_t>(c.uri, chroms, c.bin_size, c.force),
                         format, c.batch_size, c.validate_pixels);
}

static void ingest_pairs_unsorted(const LoadConfig& c) {
  assert(!c.assume_sorted);
  auto chroms = Reference::from_chrom_sizes(c.path_to_chrom_sizes);
  const auto format = format_from_string(c.format);

  const auto tmp_cooler_path = c.uri + ".tmp";

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
        {
          auto tmp_clr =
              cooler::SingleCellFile::create(tmp_cooler_path, chroms, c.bin_size, c.force);

          for (std::size_t i = 0; true; ++i) {
            SPDLOG_INFO(FMT_STRING("writing chunk #{} to intermediate file \"{}\"..."), i + 1,
                        tmp_cooler_path);
            const auto nnz = ingest_pairs_unsorted(tmp_clr.create_cell<N>(fmt::to_string(i)),
                                                   buffer, c.batch_size, format, c.validate_pixels);

            SPDLOG_INFO(FMT_STRING("done writing chunk #{} to tmp file \"{}\"."), i + 1,
                        tmp_cooler_path);
            if (nnz == 0) {
              break;
            }
          }
        }

        const cooler::SingleCellFile tmp_clr(tmp_cooler_path);
        SPDLOG_INFO(FMT_STRING("merging {} chunks into \"{}\"..."), tmp_clr.cells().size(), c.uri);
        tmp_clr.aggregate<N>(c.uri, c.force);
      },
      write_buffer);

  std::filesystem::remove(tmp_cooler_path);
}

int load_subcmd(const LoadConfig& c) {
  const auto format = format_from_string(c.format);
  const auto pixel_has_count = format == Format::COO || format == Format::BG2;

  if (c.assume_sorted && pixel_has_count) {
    SPDLOG_INFO(FMT_STRING("begin loading presorted pixels..."));
    ingest_pixels_sorted(c);
  } else if (!c.assume_sorted && pixel_has_count) {
    SPDLOG_INFO(FMT_STRING("begin loading un-sorted pixels..."));
    ingest_pixels_unsorted(c);
  } else if (c.assume_sorted && !pixel_has_count) {
    SPDLOG_INFO(FMT_STRING("begin loading presorted pairs..."));
    ingest_pairs_sorted(c);
  } else {
    SPDLOG_INFO(FMT_STRING("begin loading un-sorted pairs..."));
    ingest_pairs_unsorted(c);
  }
  return 0;
}

}  // namespace hictk::tools
