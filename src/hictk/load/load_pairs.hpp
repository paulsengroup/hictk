// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
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

  const BinTable& _bins{};  // NOLINT
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

  inline void insert_or_update(const ThinPixel<N>& pixel) {
    auto node = _buffer.find(pixel);
    if (node != _buffer.end()) {
      node->count += pixel.count;
    } else {
      _buffer.emplace(pixel);
    }
  }

  inline void read_next_batch(std::size_t batch_size) {
    _buffer.clear();

    while (!!_last_pixel) {
      const auto pixel = aggregate_pixel();
      if (!pixel) {
        break;
      }

      insert_or_update(pixel);

      if (_buffer.size() == batch_size - 1) {
        insert_or_update(_last_pixel);
        break;
      }
    }
  }
};

template <typename N>
[[nodiscard]] inline Stats ingest_pairs(
    cooler::File&& clr,  // NOLINT(*-rvalue-reference-param-not-moved)
    std::vector<ThinPixel<N>>& buffer, std::size_t batch_size, Format format, std::int64_t offset,
    bool validate_pixels) {
  buffer.reserve(batch_size);
  PairsAggregator<N>{clr.bins(), format, offset}.read_next_chunk(buffer);

  if (buffer.empty()) {
    assert(std::cin.eof());
    return {N{}, 0};
  }

  clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
  buffer.clear();

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
    std::vector<ThinPixel<float>>& buffer, Format format, std::int64_t offset) {
  const auto resolution = hf.resolutions().front();
  assert(buffer.capacity() != 0);
  buffer.reserve(buffer.capacity());
  std::size_t i = 0;

  try {
    auto t0 = std::chrono::steady_clock::now();
    for (; !std::cin.eof(); ++i) {
      PairsAggregator<float>{hf.bins(resolution), format, offset}.read_next_chunk(buffer);

      if (buffer.empty()) {
        assert(std::cin.eof());
        break;
      }
      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      t0 = t1;

      SPDLOG_INFO(FMT_STRING("preprocessing chunk #{} at {:.0f} pixels/s..."), i + 1,
                  double(buffer.size()) / delta);
      hf.add_pixels(resolution, buffer.begin(), buffer.end());
      buffer.clear();
    }
    buffer.shrink_to_fit();

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
