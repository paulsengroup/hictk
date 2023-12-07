// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

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

  const BinTable& _bins{};
  Format _format{};
  ThinPixel<N> _last_pixel{};
  std::string _line_buffer{};
  std::int64_t _offset{};

 public:
  PairsAggregator() = delete;
  inline PairsAggregator(const BinTable& bins, Format format, std::int64_t offset)
      : _bins(bins), _format(format), _offset(offset) {
    while (std::getline(std::cin, _line_buffer)) {
      if (!line_is_header(_line_buffer)) {
        break;
      }
    }
    if (!_line_buffer.empty()) {
      _last_pixel = parse_pixel<N>(_bins, _line_buffer, _format, _offset);
    }
  }

  inline bool read_next_chunk(std::vector<ThinPixel<N>>& buffer) {
    buffer.clear();
    read_next_batch(buffer.capacity());
    std::copy(_buffer.begin(), _buffer.end(), std::back_inserter(buffer));
    _buffer.clear();

    return buffer.size() == buffer.capacity();
  }

 private:
  inline ThinPixel<N> aggregate_pixel() {
    assert(!!_last_pixel);

    while (std::getline(std::cin, _line_buffer)) {
      if (_line_buffer.empty()) {
        continue;
      }

      auto p = parse_pixel<N>(_bins, _line_buffer, _format, _offset);
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

  inline void read_next_batch(std::size_t batch_size) {
    _buffer.clear();

    while (!!_last_pixel) {
      const auto pixel = aggregate_pixel();
      if (!pixel) {
        break;
      }
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
  }
};

template <typename N>
[[nodiscard]] inline std::uint64_t ingest_pairs(cooler::File&& clr,
                                                std::vector<ThinPixel<N>>& buffer,
                                                std::size_t batch_size, Format format,
                                                std::int64_t offset, bool validate_pixels) {
  buffer.reserve(batch_size);
  PairsAggregator<N>{clr.bins(), format, offset}.read_next_chunk(buffer);

  if (buffer.empty()) {
    assert(std::cin.eof());
    return {};
  }

  clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
  buffer.clear();

  clr.flush();
  return clr.nnz();
}

}  // namespace hictk::tools
