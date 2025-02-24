// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <future>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>

#include "./common.hpp"
#include "./cool.hpp"
#include "./hic.hpp"
#include "./pixel_parser.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

template <bool transpose_pixels = true, typename N>
inline void parse_pixels(PixelParser& parser, std::int64_t offset, PixelQueue<N>& queue,
                         std::atomic<bool>& early_return) {
  ThinPixel<N> buffer{};
  while (!early_return && parser.next_pixel(buffer, offset)) {
    assert(buffer.bin1_id != ThinPixel<N>::null_id);
    assert(buffer.bin2_id != ThinPixel<N>::null_id);
    assert(buffer.count != 0);

    if constexpr (transpose_pixels) {
      if (buffer.bin1_id > buffer.bin2_id) {
        std::swap(buffer.bin1_id, buffer.bin2_id);
      }
    }

    while (!queue.try_enqueue(buffer) && !early_return) {
      // NOLINTNEXTLINE(*-avoid-magic-numbers)
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  }

  if (!early_return) {
    buffer.bin1_id = ThinPixel<N>::null_id;
    buffer.bin2_id = ThinPixel<N>::null_id;
    buffer.count = 0;
    while (!queue.try_enqueue(buffer) && !early_return) {
      // NOLINTNEXTLINE(*-avoid-magic-numbers)
      std::this_thread::sleep_for(std::chrono::milliseconds(10));
    }
  }
}

template <typename N>
[[nodiscard]] inline Stats ingest_pixels(const LoadConfig& c, const BinTable& bins,
                                         std::string_view assembly, PixelQueue<N>& queue,
                                         const std::atomic<bool>& early_return) {
  if (c.output_format == "hic") {
    if constexpr (std::is_same_v<N, float>) {
      assert(c.threads > 1);
      const internal::TmpDir tmpdir{c.tmp_dir, true};
      return ingest_pixels_hic(queue, early_return, c.output_path, tmpdir(), bins.chromosomes(),
                               bins.resolution(), std::string{assembly}, c.skip_all_vs_all_matrix,
                               c.threads - 1, c.batch_size, c.compression_lvl, c.force);
    } else {
      throw std::logic_error(
          "ingest_pixels() was called with a pixel count type different than float. This is not "
          "supported when output format is .hic!");
    }
  }

  if constexpr (!std::is_same_v<N, std::int32_t> && !std::is_same_v<N, double>) {
    throw std::logic_error(
        "ingest_pixels() was called with a pixel count type different than std::int32_t and "
        "double. This is not supported when the output format is .cool!");
  } else {
    assert(c.output_format == "cool");
    const internal::TmpDir tmpdir{c.tmp_dir, true};
    const auto tmp_cooler_path =
        (tmpdir() / (std::filesystem::path{c.output_path}.filename().string() + ".tmp")).string();

    return c.assume_sorted
               ? ingest_pixels_sorted_cooler(queue, early_return, c.output_path, bins, assembly,
                                             c.batch_size, c.compression_lvl, c.force,
                                             c.validate_pixels)
               : ingest_pixels_unsorted_cooler(queue, early_return, c.output_path, tmp_cooler_path,
                                               bins, assembly, c.batch_size, c.compression_lvl,
                                               c.force, c.validate_pixels);
  }
}

template <typename N>
[[nodiscard]] inline Stats ingest_pairs(const LoadConfig& c, const BinTable& bins,
                                        std::string_view assembly, PixelQueue<N>& queue,
                                        const std::atomic<bool>& early_return) {
  if (c.output_format == "hic") {
    if constexpr (std::is_same_v<N, float>) {
      assert(c.threads > 1);
      const internal::TmpDir tmpdir{c.tmp_dir, true};
      return ingest_pairs_hic(queue, early_return, c.output_path, tmpdir(), bins.chromosomes(),
                              bins.resolution(), std::string{assembly}, c.skip_all_vs_all_matrix,
                              c.threads - 1, c.batch_size, c.compression_lvl, c.force);
    } else {
      throw std::logic_error(
          "ingest_pairs() was called with a pixel count type different than float. This is not "
          "supported when output format is .hic!");
    }
  }

  if constexpr (!std::is_same_v<N, std::int32_t> && !std::is_same_v<N, double>) {
    throw std::logic_error(
        "ingest_pixels() was called with a pixel count type different than std::int32_t and "
        "double. This is not supported when the output format is .cool!");
  } else {
    const internal::TmpDir tmpdir{c.tmp_dir, true};
    const auto tmp_cooler_path =
        (tmpdir() / (std::filesystem::path{c.output_path}.filename().string() + ".tmp")).string();

    return ingest_pairs_cooler(queue, early_return, c.output_path, tmp_cooler_path, bins, assembly,
                               c.batch_size, c.compression_lvl, c.force, c.validate_pixels);
  }
}

template <typename N>
[[nodiscard]] inline std::future<void> spawn_producer(BS::light_thread_pool& tpool,
                                                      PixelParser& parser, PixelQueue<N>& queue,
                                                      std::int64_t offset,
                                                      std::atomic<bool>& early_return,
                                                      bool transpose_lower_triangular_pixels) {
  return tpool.submit_task([&parser, &queue, offset, &early_return,
                            transpose_lower_triangular_pixels]() {
    try {
      if (transpose_lower_triangular_pixels) {
        return parse_pixels<true>(parser, offset, queue, early_return);
      }
      return parse_pixels<false>(parser, offset, queue, early_return);
    } catch (...) {
      SPDLOG_WARN(
          FMT_STRING("exception caught in thread parsing interactions: returning immediately!"));
      early_return = true;
      throw;
    }
  });
}

template <typename N>
[[nodiscard]] inline std::future<Stats> spawn_consumer(BS::light_thread_pool& tpool,
                                                       const LoadConfig& c, const BinTable& bins,
                                                       std::string_view assembly, Format format,
                                                       PixelQueue<N>& queue,
                                                       std::atomic<bool>& early_return) {
  const auto pixel_has_count = format == Format::COO || format == Format::BG2;

  return tpool.submit_task(
      [&c, &bins, assembly, &queue, &early_return, pixel_has_count]() -> Stats {
        try {
          return pixel_has_count ? ingest_pixels(c, bins, assembly, queue, early_return)
                                 : ingest_pairs(c, bins, assembly, queue, early_return);
        } catch (...) {
          SPDLOG_WARN(FMT_STRING("exception caught in thread writing interactions to file "
                                 "\"{}\": returning immediately!"),
                      c.output_path);
          early_return = true;
          throw;
        }
      });
}

}  // namespace hictk::tools
