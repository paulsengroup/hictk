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
#include "hictk/cooler.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

template <typename N>
struct PixelCmp {
  [[nodiscard]] inline bool operator()(const Pixel<N>& p1, const Pixel<N>& p2) const noexcept {
    return p1.coords < p2.coords;
  }
};

template <typename N>
class PairsAggregator {
  phmap::btree_set<Pixel<N>, PixelCmp<N>> _buffer{};
  typename phmap::btree_set<Pixel<N>, PixelCmp<N>>::const_iterator _it{};

  const BinTable& _bins{};
  Format _format{};
  Pixel<N> _last{};

 public:
  PairsAggregator() = delete;
  inline PairsAggregator(const BinTable& bins, Format format)
      : _it(_buffer.end()), _bins(bins), _format(format) {}

  inline void read_next_chunk(std::vector<Pixel<N>>& buffer) {
    buffer.clear();
    read_next_batch(buffer.capacity());
    std::copy(_buffer.begin(), _buffer.end(), std::back_inserter(buffer));
    _buffer.clear();
  }

  [[nodiscard]] inline Pixel<N> next() {
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
    auto last_bin1 = _last.coords.bin1;
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

      if (!!last_bin1 && pixel.coords.bin1 != last_bin1) {
        break;
      }
      last_bin1 = pixel.coords.bin1;
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
                                          std::vector<Pixel<N>>& buffer, std::size_t batch_size) {
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
  std::vector<Pixel<N>> buffer{};
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
inline std::string ingest_pairs_unsorted(cooler::File&& clr, std::vector<Pixel<N>>& buffer,
                                         std::size_t batch_size, Format format,
                                         bool validate_pixels) {
  buffer.reserve(batch_size);
  PairsAggregator<N>{clr.bins(), format}.read_next_chunk(buffer);

  if (buffer.empty()) {
    assert(std::cin.eof());
    return "";
  }

  clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
  buffer.clear();

  return clr.uri();
}

}  // namespace hictk::tools
