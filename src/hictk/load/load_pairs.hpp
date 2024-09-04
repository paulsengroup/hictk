// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iostream>
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
struct PixelCmp {
  [[nodiscard]] inline bool operator()(const ThinPixel<N>& p1,
                                       const ThinPixel<N>& p2) const noexcept {
    if (p1.bin1_id != p2.bin1_id) {
      return p1.bin1_id < p2.bin1_id;
    }
    return p1.bin2_id < p2.bin2_id;
  }
};

template <typename N>
class PairsAggregator {
  phmap::btree_set<ThinPixel<N>, PixelCmp<N>> _buffer{};

  PixelQueue<N>* _queue{};
  const std::atomic<bool>* _early_return{};
  ThinPixel<N> _last_pixel{};

 public:
  PairsAggregator() = delete;
  inline PairsAggregator(PixelQueue<N>& queue, const std::atomic<bool>& early_return)
      : _queue(&queue), _early_return(&early_return) {}

  inline bool read_next_chunk(std::vector<ThinPixel<N>>& buffer) {
    assert(buffer.capacity() != 0);
    buffer.clear();
    read_next_batch(buffer.capacity());
    std::copy(_buffer.begin(), _buffer.end(), std::back_inserter(buffer));
    _buffer.clear();

    return buffer.size() == buffer.capacity();
  }

 private:
  [[nodiscard]] inline ThinPixel<N> dequeue_pixel() {
    ThinPixel<N> buff{};
    while (!(*_early_return) && !_queue->wait_dequeue_timed(buff, std::chrono::milliseconds(10)));
    if (*_early_return) {
      return {ThinPixel<N>::null_id, ThinPixel<N>::null_id, 0};
    }
    return buff;
  }

  inline ThinPixel<N> aggregate_pixel() {
    while (!(*_early_return)) {
      auto p = dequeue_pixel();
      if (!p) {
        break;
      }
      if (p.bin1_id != _last_pixel.bin1_id || p.bin2_id != _last_pixel.bin2_id) {
        std::swap(p, _last_pixel);
        return p;
      }
      _last_pixel.count += p.count;
    }

    ThinPixel<N> p{};
    std::swap(p, _last_pixel);
    return p;
  }

  inline void insert_or_update(const ThinPixel<N>& pixel) {
    auto node = _buffer.find(pixel);
    if (node != _buffer.end()) {
      node->count += pixel.count;
    } else {
      _buffer.emplace(pixel);
    }
  }

  inline void read_next_batch(std::size_t batch_size) {
    assert(batch_size != 0);
    _buffer.clear();

    _last_pixel = dequeue_pixel();
    while (!!_last_pixel && _buffer.size() != batch_size - 1) {
      const auto pixel = aggregate_pixel();
      if (!pixel) {
        return;
      }

      insert_or_update(pixel);
    }

    while (!!_last_pixel && _buffer.contains(_last_pixel)) {
      insert_or_update(_last_pixel);
      _last_pixel = dequeue_pixel();
    }
    if (!!_last_pixel) {
      insert_or_update(_last_pixel);
    }
  }
};

template <typename N>
[[nodiscard]] inline Stats ingest_pairs(
    cooler::File&& clr,  // NOLINT(*-rvalue-reference-param-not-moved)
    PixelQueue<N>& queue, const std::atomic<bool>& early_return, std::vector<ThinPixel<N>>& buffer,
    std::size_t batch_size, bool validate_pixels) {
  assert(batch_size != 0);
  buffer.clear();
  buffer.reserve(batch_size);
  PairsAggregator{queue, early_return}.read_next_chunk(buffer);

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
  return {std::get<std::int64_t>(sum), nnz};
}

[[nodiscard]] inline Stats ingest_pairs(
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
      PairsAggregator{queue, early_return}.read_next_chunk(buffer);

      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      t0 = t1;

      SPDLOG_INFO(FMT_STRING("preprocessing chunk #{} at {:.0f} pixels/s..."), i + 1,
                  double(buffer.size()) / delta);
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

}  // namespace hictk::tools
