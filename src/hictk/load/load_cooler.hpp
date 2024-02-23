// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string_view>
#include <variant>
#include <vector>

#include "./common.hpp"
#include "./load_pairs.hpp"
#include "./load_pixels.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

inline Stats ingest_pixels_unsorted_cooler(std::string_view uri, std::string_view tmp_cooler_path,
                                           const Reference& chromosomes, std::uint32_t bin_size,
                                           std::string_view assembly, std::int64_t offset,
                                           Format format, std::size_t batch_size,
                                           std::uint32_t compression_lvl, bool force,
                                           bool count_as_float, bool validate_pixels) {
  SPDLOG_INFO(FMT_STRING("begin loading unsorted pixels into a .cool file..."));
  const BinTable bins(chromosomes, bin_size);
  PixelBuffer write_buffer{};
  if (count_as_float) {
    write_buffer = FPBuff(batch_size);
  } else {
    write_buffer = IntBuff(batch_size);
  }

  const auto stats = std::visit(
      [&](auto& buffer) {
        using N = decltype(buffer.front().count);
        Stats local_stats{N{}, 0};
        {
          auto sclr_attrs = cooler::SingleCellAttributes::init(bins.bin_size());
          sclr_attrs.assembly = assembly;

          auto tmp_clr = cooler::SingleCellFile::create(tmp_cooler_path, bins, force, sclr_attrs);
          auto attrs = cooler::Attributes::init(bins.bin_size());
          attrs.assembly = assembly;

          for (std::size_t i = 0; true; ++i) {
            SPDLOG_INFO(FMT_STRING("writing chunk #{} to intermediate file \"{}\"..."), i + 1,
                        tmp_cooler_path);
            const auto partial_stats = ingest_pixels_unsorted(
                tmp_clr.create_cell<N>(fmt::to_string(i), attrs,
                                       cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl),
                buffer, format, offset, validate_pixels);
            local_stats += partial_stats;
            SPDLOG_INFO(FMT_STRING("done writing chunk #{} to tmp file \"{}\"."), i + 1,
                        tmp_cooler_path);
            if (partial_stats.nnz == 0) {
              break;
            }
          }
        }
        const cooler::SingleCellFile tmp_clr(tmp_cooler_path);
        SPDLOG_INFO(FMT_STRING("merging {} chunks into \"{}\"..."), tmp_clr.cells().size(), uri);
        tmp_clr.aggregate<N>(uri, force, compression_lvl);

        return local_stats;
      },
      write_buffer);
  std::filesystem::remove(tmp_cooler_path);

  return stats;
}

inline Stats ingest_pixels_sorted_cooler(std::string_view uri, const Reference& chromosomes,
                                         std::uint32_t bin_size, std::string_view assembly,
                                         std::int64_t offset, Format format, std::size_t batch_size,
                                         std::uint32_t compression_lvl, bool force,
                                         bool count_as_float, bool validate_pixels) {
  SPDLOG_INFO(FMT_STRING("begin loading pre-sorted pixels into a .cool file..."));
  auto attrs = cooler::Attributes::init(bin_size);
  attrs.assembly = assembly;
  if (count_as_float) {
    return ingest_pixels_sorted<double>(
        cooler::File::create<double>(uri, chromosomes, bin_size, force, attrs,
                                     cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl),
        format, offset, batch_size, validate_pixels);
  }
  return ingest_pixels_sorted<std::int32_t>(
      cooler::File::create<std::int32_t>(uri, chromosomes, bin_size, force, attrs,
                                         cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl),
      format, offset, batch_size, validate_pixels);
}

inline Stats ingest_pairs_cooler(std::string_view uri, std::string_view tmp_cooler_path,
                                 const BinTable& bins, std::string_view assembly,
                                 std::int64_t offset, Format format, std::size_t batch_size,
                                 std::uint32_t compression_lvl, bool force, bool count_as_float,
                                 bool validate_pixels) {
  PixelBuffer write_buffer{};
  if (count_as_float) {
    write_buffer = FPBuff{};
  } else {
    write_buffer = IntBuff{};
  }

  std::visit(
      [&](auto& buffer) {
        using N = decltype(buffer.begin()->count);
        {
          auto sclr_attrs = cooler::SingleCellAttributes::init(bins.bin_size());
          sclr_attrs.assembly = assembly;

          auto tmp_clr = cooler::SingleCellFile::create(tmp_cooler_path, bins, force, sclr_attrs);
          auto attrs = cooler::Attributes::init(bins.bin_size());
          attrs.assembly = assembly;

          for (std::size_t i = 0; true; ++i) {
            SPDLOG_INFO(FMT_STRING("writing chunk #{} to intermediate file \"{}\"..."), i + 1,
                        tmp_cooler_path);
            const auto partial_stats = ingest_pairs(
                tmp_clr.create_cell<N>(fmt::to_string(i), attrs,
                                       cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl),
                buffer, batch_size, format, offset, validate_pixels);

            SPDLOG_INFO(FMT_STRING("done writing chunk #{} to tmp file \"{}\"."), i + 1,
                        tmp_cooler_path);
            if (partial_stats.nnz == 0) {
              break;
            }
          }
        }

        const cooler::SingleCellFile tmp_clr(tmp_cooler_path);
        SPDLOG_INFO(FMT_STRING("merging {} chunks into \"{}\"..."), tmp_clr.cells().size(), uri);
        tmp_clr.aggregate<N>(uri, force);
      },
      write_buffer);

  std::filesystem::remove(tmp_cooler_path);

  const cooler::File clr(uri);
  const auto nnz = clr.nnz();
  const auto sum = clr.attributes().sum.value();

  if (clr.has_float_pixels()) {
    return {std::get<double>(sum), nnz};
  }
  return {std::get<std::int64_t>(sum), nnz};
}

}  // namespace hictk::tools
