// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>

#include "./load.hpp"
#include "./pairs_aggregator.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk::tools {

[[nodiscard]] static Stats ingest_pairs(
    hic::internal::HiCFileWriter&& hf,  // NOLINT(*-rvalue-reference-param-not-moved)
    PixelQueue<float>& queue, const std::atomic<bool>& early_return,
    std::vector<ThinPixel<float>>& buffer, std::size_t batch_size) {
  const auto resolution = hf.resolutions().front();
  assert(batch_size != 0);
  buffer.clear();
  buffer.reserve(batch_size);
  std::size_t i = 0;

  try {
    auto t0 = std::chrono::steady_clock::now();
    for (; !early_return; ++i) {
      buffer.clear();
      // clang-8 does not like when read_next_chunk() is called directly when PairsAggregator is
      // constructed
      PairsAggregator aggr(queue, early_return);
      aggr.read_next_chunk(buffer);

      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      t0 = t1;

      SPDLOG_INFO(FMT_STRING("preprocessing chunk #{} at {:.0f} pixels/s..."), i + 1,
                  static_cast<double>(buffer.size()) / delta);
      hf.add_pixels(resolution, buffer.begin(), buffer.end());

      if (buffer.size() != buffer.capacity()) {
        break;
      }
    }

    hf.serialize();
    const auto stats = hf.stats(resolution);
    return {stats.sum, stats.nnz};
  } catch (const std::exception& e) {
    const auto i0 = i * buffer.capacity();
    const auto i1 = i0 + buffer.size();
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while processing chunk {}-{}: {}"), i0, i1, e.what()));
  }
}

[[nodiscard]] static Stats ingest_pixels(
    hic::internal::HiCFileWriter&& hf,  // NOLINT(*-rvalue-reference-param-not-moved)
    PixelQueue<float>& queue, const std::atomic<bool>& early_return,
    std::vector<ThinPixel<float>>& buffer) {
  assert(buffer.capacity() != 0);

  std::size_t i = 0;
  Stats stats{0.0, 0};
  try {
    auto t0 = std::chrono::steady_clock::now();
    const auto& bins = hf.bins(hf.resolutions().front());
    for (; !early_return; ++i) {
      stats += read_batch(queue, early_return, buffer);

      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      t0 = t1;
      SPDLOG_INFO(FMT_STRING("preprocessing chunk #{} at {:.0f} pixels/s..."), i + 1,
                  static_cast<double>(buffer.size()) / delta);
      hf.add_pixels(bins.resolution(), buffer.begin(), buffer.end());
      if (buffer.size() != buffer.capacity()) {
        break;
      }
      buffer.clear();
    }
    hf.serialize();
    return stats;
  } catch (const std::exception& e) {
    const auto i0 = i * buffer.capacity();
    const auto i1 = i0 + buffer.size();
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while processing chunk {}-{}: {}"), i0, i1, e.what()));
  }
}

Stats ingest_pixels_hic(PixelQueue<float>& pixel_queue, const std::atomic<bool>& early_return,
                        std::string_view uri, const std::filesystem::path& tmp_dir,
                        const Reference& chromosomes, std::uint32_t bin_size,
                        const std::string& assembly, bool skip_all_vs_all_matrix,
                        std::size_t threads, std::size_t batch_size, std::uint32_t compression_lvl,
                        bool force) {
  SPDLOG_INFO("begin loading pixels into a .hic file...");

  if (force) {
    std::filesystem::remove(uri);  // NOLINT
  }

  hic::internal::HiCFileWriter hf(uri, chromosomes, {bin_size}, assembly, threads, batch_size,
                                  tmp_dir, compression_lvl, skip_all_vs_all_matrix);

  std::vector<ThinPixel<float>> write_buffer(batch_size);
  return ingest_pixels(std::move(hf), pixel_queue, early_return, write_buffer);
}

Stats ingest_pairs_hic(PixelQueue<float>& pixel_queue, const std::atomic<bool>& early_return,
                       std::string_view uri, const std::filesystem::path& tmp_dir,
                       const Reference& chromosomes, std::uint32_t bin_size,
                       const std::string& assembly, bool skip_all_vs_all_matrix,
                       std::size_t threads, std::size_t batch_size, std::uint32_t compression_lvl,
                       bool force) {
  if (force) {
    std::filesystem::remove(uri);  // NOLINT
  }

  hic::internal::HiCFileWriter hf(uri, chromosomes, {bin_size}, assembly, threads, batch_size,
                                  tmp_dir, compression_lvl, skip_all_vs_all_matrix);

  std::vector<ThinPixel<float>> buffer(batch_size);
  return ingest_pairs(std::move(hf), pixel_queue, early_return, buffer, batch_size);
}

}  // namespace hictk::tools
