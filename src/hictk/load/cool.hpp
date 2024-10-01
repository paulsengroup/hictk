// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#include "./common.hpp"
#include "./pairs_aggregator.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/singlecell_cooler.hpp"
#include "hictk/pixel.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

template <typename N>
[[nodiscard]] inline Stats ingest_pairs(
    cooler::File&& clr,  // NOLINT(*-rvalue-reference-param-not-moved)
    PixelQueue<N>& queue, const std::atomic<bool>& early_return, std::vector<ThinPixel<N>>& buffer,
    std::size_t batch_size, bool validate_pixels) {
  assert(batch_size != 0);
  buffer.clear();
  buffer.reserve(batch_size);
  // clang-8 does not like when read_next_chunk() is called directly when PairsAggregator is
  // constructed
  PairsAggregator aggr(queue, early_return);
  aggr.read_next_chunk(buffer);

  if (buffer.empty()) {
    return {N{}, 0};
  }

  clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);

  clr.flush();
  const auto nnz = clr.nnz();
  const auto sum = clr.attributes().sum.value();

  if (clr.has_float_pixels()) {
    return {std::get<double>(sum), nnz};
  }
  return {static_cast<std::uint64_t>(std::get<std::int64_t>(sum)), nnz};
}

template <typename N>
[[nodiscard]] inline Stats ingest_pixels_sorted(
    cooler::File&& clr,  // NOLINT(*-rvalue-reference-param-not-moved)
    PixelQueue<N>& queue, const std::atomic<bool>& early_return, std::size_t batch_size,
    bool validate_pixels) {
  std::vector<ThinPixel<N>> buffer(batch_size);

  std::size_t i = 0;
  Stats stats{N{}, 0};
  try {
    for (; true; ++i) {
      SPDLOG_INFO(FMT_STRING("processing chunk #{}..."), i + 1);
      stats += read_batch(queue, early_return, buffer);

      clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
      if (buffer.size() != batch_size) {
        return stats;
      }
      buffer.clear();
    }
  } catch (const std::exception& e) {
    const auto i0 = i * buffer.capacity();
    const auto i1 = i0 + buffer.size();
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while processing chunk {}-{}: {}"), i0, i1, e.what()));
  }
}

template <typename N>
[[nodiscard]] inline Stats ingest_pixels_unsorted(
    cooler::File&& clr,  // NOLINT(*-rvalue-reference-param-not-moved)
    PixelQueue<N>& queue, const std::atomic<bool>& early_return, std::vector<ThinPixel<N>>& buffer,
    bool validate_pixels) {
  assert(buffer.capacity() != 0);

  auto stats = read_batch(queue, early_return, buffer);

  if (buffer.empty()) {
    return {N{}, 0};
  }

  std::sort(buffer.begin(), buffer.end());
  clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);

  clr.flush();
  return stats;
}

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

  {
    const cooler::SingleCellFile tmp_clr(tmp_cooler_path);
    SPDLOG_INFO(FMT_STRING("merging {} chunks into \"{}\"..."), tmp_clr.cells().size(), uri);
    std::ignore = tmp_clr.aggregate<N>(uri, force, compression_lvl);
  }

  std::filesystem::remove(tmp_cooler_path);  // NOLINT

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

  {
    const cooler::SingleCellFile tmp_clr(tmp_cooler_path);
    SPDLOG_INFO(FMT_STRING("merging {} chunks into \"{}\"..."), tmp_clr.cells().size(), uri);
    std::ignore = tmp_clr.aggregate<N>(uri, force);
  }

  std::filesystem::remove(tmp_cooler_path);  // NOLINT

  const cooler::File clr(uri);
  const auto nnz = clr.nnz();
  const auto sum = clr.attributes().sum.value();

  if (clr.has_float_pixels()) {
    return {std::get<double>(sum), nnz};
  }
  return {static_cast<std::uint64_t>(std::get<std::int64_t>(sum)), nnz};
}

template <typename N>
static Stats ingest_pixels_cooler(const LoadConfig& c, const BinTable& bins,
                                  std::string_view assembly, PixelQueue<N>& queue,
                                  const std::atomic<bool>& early_return) {
  assert(c.output_format == "cool");
  const internal::TmpDir tmpdir{c.tmp_dir, true};
  const auto tmp_cooler_path =
      (tmpdir() / (std::filesystem::path{c.output_path}.filename().string() + ".tmp")).string();

  return c.assume_sorted ? ingest_pixels_sorted_cooler(queue, early_return, c.output_path, bins,
                                                       assembly, c.batch_size, c.compression_lvl,
                                                       c.force, c.validate_pixels)
                         : ingest_pixels_unsorted_cooler(
                               queue, early_return, c.output_path, tmp_cooler_path, bins, assembly,
                               c.batch_size, c.compression_lvl, c.force, c.validate_pixels);
}

}  // namespace hictk::tools
