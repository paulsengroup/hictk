// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <chrono>
#include <cstddef>
#include <vector>

#include "./common.hpp"
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

}  // namespace hictk::tools
