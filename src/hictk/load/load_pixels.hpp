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
#include <exception>
#include <stdexcept>
#include <string>
#include <vector>

#include "./common.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

template <typename N>
inline Stats read_batch(PixelQueue<N>& queue, const std::atomic<bool>& early_return,
                        std::vector<ThinPixel<N>>& buffer) {
  assert(buffer.capacity() != 0);
  buffer.clear();
  Stats stats{N{}, 0};
  ThinPixel<N> pixel{};

  while (!early_return) {
    if (!queue.wait_dequeue_timed(pixel, std::chrono::milliseconds(10))) {
      continue;
    }

    if (pixel.bin1_id == ThinPixel<N>::null_id && pixel.bin2_id == ThinPixel<N>::null_id &&
        pixel.count == 0) {
      // EOQ signal received
      return stats;
    }

    buffer.emplace_back(pixel);
    stats.nnz++;
    if constexpr (std::is_floating_point_v<N>) {
      std::get<double>(stats.sum) += conditional_static_cast<double>(pixel.count);
    } else {
      std::get<std::uint64_t>(stats.sum) += conditional_static_cast<std::uint64_t>(pixel.count);
    }
    if (buffer.size() == buffer.capacity()) {
      return stats;
    }
  }

  return stats;
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

[[nodiscard]] inline Stats ingest_pixels(
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
                  double(buffer.size()) / delta);
      hf.add_pixels(bins.resolution(), buffer.begin(), buffer.end());
      if (buffer.size() != buffer.capacity()) {
        break;
      }
      buffer.clear();
    }
    hf.serialize();
    assert(buffer.empty());
    return stats;
  } catch (const std::exception& e) {
    const auto i0 = i * buffer.capacity();
    const auto i1 = i0 + buffer.size();
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while processing chunk {}-{}: {}"), i0, i1, e.what()));
  }
}

}  // namespace hictk::tools
