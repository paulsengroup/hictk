// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename PixelIt>
inline CoarsenPixels<PixelIt>::CoarsenPixels(PixelIt first_pixel, PixelIt last_pixel,
                                             std::shared_ptr<const BinTable> source_bins,
                                             std::size_t factor)
    : _first(std::move(first_pixel)),
      _last(std::move(last_pixel)),
      _src_bins(std::move(source_bins)),
      _dest_bins(std::make_shared<const BinTable>(_src_bins->chromosomes(),
                                                  _src_bins->bin_size() * factor)),
      _factor(factor) {
  if (factor < 2) {
    throw std::logic_error("coarsening factor should be > 1");
  }
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::begin() const -> iterator {
  return iterator{_first, _last, _src_bins, _dest_bins};
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::cbegin() const -> iterator {
  return this->begin();
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::end() const -> iterator {
  return iterator::at_end(_last, _src_bins, _dest_bins);
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::cend() const -> iterator {
  return this->end();
}

template <typename PixelIt>
inline const BinTable &CoarsenPixels<PixelIt>::src_bins() const noexcept {
  return *this->src_bins_ptr();
}
template <typename PixelIt>
inline const BinTable &CoarsenPixels<PixelIt>::dest_bins() const noexcept {
  return *this->dest_bins_ptr();
}
template <typename PixelIt>
inline std::shared_ptr<const BinTable> CoarsenPixels<PixelIt>::src_bins_ptr() const noexcept {
  return this->_src_bins;
}
template <typename PixelIt>
inline std::shared_ptr<const BinTable> CoarsenPixels<PixelIt>::dest_bins_ptr() const noexcept {
  return this->_dest_bins;
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::read_all() const -> std::vector<ThinPixel<N>> {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<ThinPixel<N>> buff{};
  std::copy(this->begin(), this->end(), std::back_inserter(buff));
  return buff;
}

template <typename PixelIt>
inline CoarsenPixels<PixelIt>::iterator::iterator(PixelIt first, PixelIt last,
                                                  std::shared_ptr<const BinTable> src_bins,
                                                  std::shared_ptr<const BinTable> dest_bins)
    : _pixel_it(std::move(first)),
      _pixel_last(std::move(last)),
      _src_bins(std::move(src_bins)),
      _dest_bins(std::move(dest_bins)),
      _buffer(std::make_shared<BufferT>()) {
  assert(_dest_bins->bin_size() > _src_bins->bin_size());

  if (_pixel_it == _pixel_last) {
    *this = at_end(_pixel_last, _src_bins, _dest_bins);
  } else {
    process_next_row();
  }
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::at_end(PixelIt last,
                                                     std::shared_ptr<const BinTable> src_bins,
                                                     std::shared_ptr<const BinTable> dest_bins)
    -> iterator {
  iterator it{};
  it._pixel_it = last;
  it._pixel_last = last;
  it._src_bins = std::move(src_bins);
  it._dest_bins = std::move(dest_bins);
  return it;
}

template <typename PixelIt>
inline void CoarsenPixels<PixelIt>::iterator::process_next_row() {
  if (_pixel_it == _pixel_last) {
    *this = at_end(_pixel_last, _src_bins, _dest_bins);
    return;
  }
  auto buffer = coarsen_chunk_pass1();
  coarsen_chunk_pass2(buffer);
  if (!!_buffer) {
    _it = _buffer->begin();
  }
}

// Loop over the current chunk and coarse pixels based on bin1_ids
template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::coarsen_chunk_pass1() -> ColumnMerger {
  assert(_pixel_it != _pixel_last);
  ColumnMerger merger{};

  // Compute the first and last bins mapping to chunk we are processing
  const auto factor = _dest_bins->bin_size() / _src_bins->bin_size();
  _bin1_id_chunk_start = (_src_bins->at(_pixel_it->bin1_id).rel_id() / factor) * factor;
  _bin1_id_chunk_end = _bin1_id_chunk_start + factor;

  // Init one empty buffer for each row in the current chunk
  for (std::size_t i = 0; i < factor; ++i) {
    merger.emplace(_bin1_id_chunk_start + i, BufferT{});
  }

  // Init current buffer
  auto current_bin = _bin1_id_chunk_start;
  auto *row_buffer = &merger.at(current_bin);

  for (std::size_t i = 0; _pixel_it != _pixel_last; ++i) {
    // We need to map pixel coordinates instead of just dividing bin ids by the coarsening factor to
    // avoid mapping the last pixel in a chromosome i and the first pixel in chromosome i+1 as to
    // the same coarse bin
    const PixelCoordinates src_coords{_src_bins->at(_pixel_it->bin1_id),
                                      _src_bins->at(_pixel_it->bin2_id)};
    const PixelCoordinates dest_coords{_dest_bins->at(src_coords.bin1.interval()).first,
                                       _dest_bins->at(src_coords.bin2.interval()).first};

    if (const auto id = src_coords.bin1.rel_id();
        id < _bin1_id_chunk_start || id >= _bin1_id_chunk_end) {
      if (i == 0) {
        return coarsen_chunk_pass1();  // found an empty row
      }
      break;  // done processing current chunk
    }

    auto pixel = ThinPixel<N>{dest_coords.bin1.id(), dest_coords.bin2.id(), _pixel_it->count};

    if (src_coords.bin1.rel_id() != current_bin) {
      // We're processing a new row: update row buffer
      current_bin = src_coords.bin1.rel_id();
      row_buffer = &merger.at(current_bin);
    }

    // Insert or merge pixel
    if (row_buffer->empty() || row_buffer->back().bin2_id != pixel.bin2_id) {
      row_buffer->emplace_back(pixel);
    } else {
      row_buffer->back().count += pixel.count;
    }
    ++_pixel_it;
  }

  return merger;
}

// Loop over pixels produced by coarsen_chunk_pass1 and merge pixel sharing the same coords
template <typename PixelIt>
inline void CoarsenPixels<PixelIt>::iterator::coarsen_chunk_pass2(const ColumnMerger &col_merger) {
  if (col_merger.empty()) {
    *this = at_end(_pixel_last, _src_bins, _dest_bins);
    return;
  }

  std::vector<RowIt> heads{};
  std::vector<RowIt> tails{};

  for (const auto &[_, buff] : col_merger) {
    if (buff.begin() != buff.end()) {
      heads.emplace_back(buff.begin());
      tails.emplace_back(buff.end());
    }
  }

  using PixelMerger = internal::PixelMerger<RowIt>;
  PixelMerger merger(heads, tails);

  if (_buffer.use_count() != 1) {
    _buffer = std::make_shared<BufferT>();
  }
  _buffer->clear();
  _buffer->emplace_back(merger.next());
  while (!!_buffer->back()) {
    const auto value = merger.next();
    if (!value) {
      break;  // all pixels have been processed
    }

    // Insert or merge pixels
    auto &last_pixel = _buffer->back();
    if (value.bin1_id == last_pixel.bin1_id && value.bin2_id == last_pixel.bin2_id) {
      last_pixel.count += value.count;
    } else {
      _buffer->emplace_back(value);
    }
  }
}

template <typename PixelIt>
inline bool CoarsenPixels<PixelIt>::iterator::operator==(
    const CoarsenPixels::iterator &other) const noexcept {
  return this->_buffer == other._buffer && this->_it == other._it &&
         this->_src_bins == other._src_bins && this->_dest_bins == other._dest_bins;
}

template <typename PixelIt>
inline bool CoarsenPixels<PixelIt>::iterator::operator!=(
    const CoarsenPixels::iterator &other) const noexcept {
  return !(*this == other);
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::operator*() const -> const_reference {
  return *this->_it;
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::operator->() const -> const_pointer {
  return &(**this);
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::operator++() -> iterator & {
  assert(_buffer);

  if (++_it == _buffer->end()) {
    this->process_next_row();
  }
  return *this;
}

template <typename PixelIt>
inline auto CoarsenPixels<PixelIt>::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

}  // namespace hictk::transformers
