// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "./common.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
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
  typename phmap::btree_set<ThinPixel<N>, PixelCmp<N>>::const_iterator _it{};

  const BinTable& _bins{};
  Format _format{};
  ThinPixel<N> _last{};

 public:
  PairsAggregator() = delete;
  inline PairsAggregator(const BinTable& bins, Format format)
      : _it(_buffer.end()), _bins(bins), _format(format) {}

  inline void read_next_chunk(std::vector<ThinPixel<N>>& buffer) {
    buffer.clear();
    read_next_batch(buffer.capacity());
    std::copy(_buffer.begin(), _buffer.end(), std::back_inserter(buffer));
    _buffer.clear();
  }

  [[nodiscard]] inline ThinPixel<N> next() {
    if (_it == _buffer.end()) {
      read_next_batch();
    }
    if (_buffer.empty()) {
      auto last = _last;
      _last = {};
      return last;
    }
    return *_it++;
  }

 private:
  inline void read_next_batch() {
    auto last_bin1 = _last.bin1_id;
    std::string line{};

    _buffer.clear();
    if (!!_last) {
      _buffer.emplace(_last);
    }

    while (std::getline(std::cin, line)) {
      if (line_is_header(line)) {
        continue;
      }

      auto pixel = parse_pixel<N>(_bins, line, _format);
      auto node = _buffer.find(pixel);
      if (node != _buffer.end()) {
        node->count += pixel.count;
      } else {
        _buffer.emplace(pixel);
      }

      if (last_bin1 != ThinPixel<N>::null_id && pixel.bin1_id != last_bin1) {
        break;
      }
      last_bin1 = pixel.bin1_id;
    }

    _it = _buffer.begin();
    if (_buffer.empty()) {
      _last = {};
    } else {
      _last = *_buffer.rbegin();
      _buffer.erase(_last);
    }
  }

  inline void read_next_batch(std::size_t batch_size) {
    std::string line{};

    _buffer.clear();
    while (std::getline(std::cin, line)) {
      if (line_is_header(line)) {
        continue;
      }

      auto pixel = parse_pixel<N>(_bins, line, _format);
      auto node = _buffer.find(pixel);
      if (node != _buffer.end()) {
        node->count += pixel.count;
      } else {
        _buffer.emplace(pixel);
      }

      if (_buffer.size() == batch_size) {
        break;
      }
    }

    _it = _buffer.begin();
  }
};

template <typename N>
inline void read_sort_and_aggregate_batch(PairsAggregator<N>& aggregator,
                                          std::vector<ThinPixel<N>>& buffer,
                                          std::size_t batch_size) {
  buffer.reserve(batch_size);
  buffer.clear();

  while (true) {
    if (buffer.size() == batch_size) {
      return;
    }

    auto pixel = aggregator.next();
    if (!pixel) {
      break;
    }
    buffer.emplace_back(std::move(pixel));
  }
}

template <typename N>
inline void ingest_pairs_sorted(cooler::File&& clr, Format format, std::size_t batch_size,
                                bool validate_pixels) {
  PairsAggregator<N> aggregator{clr.bins(), format};
  std::vector<ThinPixel<N>> buffer{};
  buffer.reserve(batch_size);

  std::size_t i0 = 0;
  std::size_t i1 = i0;
  try {
    for (; true; ++i1) {
      if (buffer.size() == batch_size) {
        spdlog::info(FMT_STRING("processing chunk {}-{}..."), i0, i1);
        clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
        buffer.clear();
        i0 = i1;
      }

      auto pixel = aggregator.next();
      if (!pixel) {
        break;
      }
      buffer.emplace_back(std::move(pixel));
    }

    if (!buffer.empty()) {
      clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while processing chunk {}-{}: {}"), i0, i1, e.what()));
  }
}

template <typename N>
[[nodiscard]] inline CoolerChunk<N> ingest_pairs_unsorted(cooler::File&& clr,
                                                          std::vector<ThinPixel<N>>& buffer,
                                                          std::size_t batch_size, Format format,
                                                          bool validate_pixels) {
  buffer.reserve(batch_size);
  PairsAggregator<N>{clr.bins(), format}.read_next_chunk(buffer);

  if (buffer.empty()) {
    assert(std::cin.eof());
    return {};
  }

  clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
  buffer.clear();

  clr.flush();
  auto sel = clr.fetch();
  return {clr.uri(), sel.begin<N>(), sel.end<N>()};
}

}  // namespace hictk::tools
