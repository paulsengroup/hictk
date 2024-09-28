// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#if __has_include(<readerwriterqueue.h>)
#include <readerwriterqueue.h>
#else
#include <readerwriterqueue/readerwriterqueue.h>
#endif

#include <atomic>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <string_view>
#include <variant>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

template <typename N>
using PixelQueue = moodycamel::BlockingReaderWriterQueue<ThinPixel<N>>;

using PixelQueueVar = std::variant<PixelQueue<std::int32_t>, PixelQueue<float>, PixelQueue<double>>;

using IntBuff = std::vector<ThinPixel<std::int32_t>>;
using FPBuff = std::vector<ThinPixel<double>>;
using PixelBuffer = std::variant<IntBuff, FPBuff>;

struct Stats {
  std::variant<std::uint64_t, double> sum{0.0};
  std::uint64_t nnz{};

  template <typename N>
  inline Stats(N sum_, std::uint64_t nnz_) : nnz(nnz_) {
    if constexpr (std::is_floating_point_v<N>) {
      sum = conditional_static_cast<double>(sum_);
    } else {
      sum = conditional_static_cast<std::uint64_t>(sum_);
    }
  }

  Stats& operator+=(const Stats& other);
};

enum class Format : std::uint_fast8_t { COO, BG2, VP, _4DN };
[[nodiscard]] Format format_from_string(std::string_view s);

template <typename N>
[[nodiscard]] inline Stats read_batch(PixelQueue<N>& queue, const std::atomic<bool>& early_return,
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

}  // namespace hictk::tools
