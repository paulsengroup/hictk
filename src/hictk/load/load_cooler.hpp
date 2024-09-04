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

template <typename N>
[[nodiscard]] inline Stats ingest_pixels_unsorted_cooler(
    PixelQueue<N>& pixel_queue, const std::atomic<bool>& early_return, std::string_view uri,
    std::string_view tmp_cooler_path, const BinTable& bins, std::string_view assembly,
    std::size_t batch_size, std::uint32_t compression_lvl, bool force, bool validate_pixels) {
  static_assert(std::is_same_v<N, std::int32_t> || std::is_same_v<N, double>);
  SPDLOG_INFO(FMT_STRING("begin loading unsorted pixels into a .cool file..."));
  Stats stats{N{}, 0};
  std::vector<ThinPixel<N>> write_buffer(batch_size);

  {
    auto sclr_attrs = cooler::SingleCellAttributes::init(bins.resolution());
    sclr_attrs.assembly = assembly;

    auto tmp_clr = cooler::SingleCellFile::create(tmp_cooler_path, bins, force, sclr_attrs);
    auto attrs = cooler::Attributes::init(bins.resolution());
    attrs.assembly = assembly;

    for (std::size_t i = 0; true; ++i) {
      SPDLOG_INFO(FMT_STRING("writing chunk #{} to intermediate file \"{}\"..."), i + 1,
                  tmp_cooler_path);
      const auto partial_stats = ingest_pixels_unsorted(
          tmp_clr.create_cell<N>(fmt::to_string(i), attrs, cooler::DEFAULT_HDF5_CACHE_SIZE * 4,
                                 compression_lvl),
          pixel_queue, early_return, write_buffer, validate_pixels);
      stats += partial_stats;
      SPDLOG_INFO(FMT_STRING("done writing chunk #{} to tmp file \"{}\"."), i + 1, tmp_cooler_path);
      if (write_buffer.size() != batch_size) {
        break;
      }
    }
  }
  const cooler::SingleCellFile tmp_clr(tmp_cooler_path);
  SPDLOG_INFO(FMT_STRING("merging {} chunks into \"{}\"..."), tmp_clr.cells().size(), uri);
  std::ignore = tmp_clr.aggregate<N>(uri, force, compression_lvl);

  std::filesystem::remove(tmp_cooler_path);

  return stats;
}

template <typename N>
[[nodiscard]] inline Stats ingest_pixels_sorted_cooler(
    PixelQueue<N>& pixel_queue, const std::atomic<bool>& early_return, std::string_view uri,
    const BinTable& bins, std::string_view assembly, std::size_t batch_size,
    std::uint32_t compression_lvl, bool force, bool validate_pixels) {
  static_assert(std::is_same_v<N, std::int32_t> || std::is_same_v<N, double>);
  SPDLOG_INFO(FMT_STRING("begin loading pre-sorted pixels into a .cool file..."));
  auto attrs = cooler::Attributes::init(bins.resolution());
  attrs.assembly = assembly;
  return ingest_pixels_sorted<N>(
      cooler::File::create<N>(uri, bins, force, attrs, cooler::DEFAULT_HDF5_CACHE_SIZE * 4,
                              compression_lvl),
      pixel_queue, early_return, batch_size, validate_pixels);
}

template <typename N>
[[nodiscard]] inline Stats ingest_pairs_cooler(
    PixelQueue<N>& pixel_queue, const std::atomic<bool>& early_return, std::string_view uri,
    std::string_view tmp_cooler_path, const BinTable& bins, std::string_view assembly,
    std::size_t batch_size, std::uint32_t compression_lvl, bool force, bool validate_pixels) {
  static_assert(std::is_same_v<N, std::int32_t> || std::is_same_v<N, double>);
  SPDLOG_INFO(FMT_STRING("begin loading pairwise interactions into a .cool file..."));
  std::vector<ThinPixel<N>> write_buffer(batch_size);
  {
    auto sclr_attrs = cooler::SingleCellAttributes::init(bins.resolution());
    sclr_attrs.assembly = assembly;

    auto tmp_clr = cooler::SingleCellFile::create(tmp_cooler_path, bins, force, sclr_attrs);
    auto attrs = cooler::Attributes::init(bins.resolution());
    attrs.assembly = assembly;

    for (std::size_t i = 0; true; ++i) {
      SPDLOG_INFO(FMT_STRING("writing chunk #{} to intermediate file \"{}\"..."), i + 1,
                  tmp_cooler_path);
      const auto partial_stats =
          ingest_pairs(tmp_clr.create_cell<N>(fmt::to_string(i), attrs,
                                              cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl),
                       pixel_queue, early_return, write_buffer, batch_size, validate_pixels);
      SPDLOG_INFO(FMT_STRING("done writing chunk #{} to tmp file \"{}\"."), i + 1, tmp_cooler_path);
      if (write_buffer.size() != batch_size || partial_stats.nnz == 0) {
        break;
      }
    }
  }

  const cooler::SingleCellFile tmp_clr(tmp_cooler_path);
  SPDLOG_INFO(FMT_STRING("merging {} chunks into \"{}\"..."), tmp_clr.cells().size(), uri);
  std::ignore = tmp_clr.aggregate<N>(uri, force);

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
