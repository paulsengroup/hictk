// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <memory>
#include <utility>

#include "hictk/numeric_utils.hpp"

namespace hictk::cooler {

inline PixelSelector::PixelSelector(std::shared_ptr<const Index> index,
                                    const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                    const Dataset &pixels_count, PixelCoordinates coords,
                                    std::shared_ptr<const balancing::Weights> weights) noexcept
    : PixelSelector(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count, coords, coords,
                    std::move(weights)) {}

inline PixelSelector::PixelSelector(std::shared_ptr<const Index> index,
                                    const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                    const Dataset &pixels_count, PixelCoordinates coord1,
                                    PixelCoordinates coord2,
                                    std::shared_ptr<const balancing::Weights> weights) noexcept
    : _coord1(std::move(coord1)),
      _coord2(std::move(coord2)),
      _index(std::move(index)),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count),
      _weights(std::move(weights)) {
  assert(_index);
}

inline PixelSelector::PixelSelector(std::shared_ptr<const Index> index,
                                    const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                    const Dataset &pixels_count,
                                    std::shared_ptr<const balancing::Weights> weights) noexcept
    : _index(std::move(index)),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count),
      _weights(std::move(weights)) {}

inline bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  // clang-format off
  return this->begin<int>() == other.begin<int>() &&
         this->end<int>() == other.end<int>() &&
         this->_weights == other._weights;
  // clang-format on
}

inline bool PixelSelector::operator!=(const PixelSelector &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelector::begin() const -> iterator<N> {
  return this->cbegin<N>();
}

template <typename N>
inline auto PixelSelector::end() const -> iterator<N> {
  return this->cend<N>();
}

template <typename N>
inline auto PixelSelector::cbegin() const -> iterator<N> {
  if constexpr (std::is_integral_v<N>) {
    if (!!_weights) {
      throw std::logic_error(
          "iterator template parameter should be of floating point type when processing balanced "
          "matrices.");
    }
  }

  if (!this->_coord1) {
    assert(!this->_coord2);
    return iterator<N>{this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                       *this->_pixels_count, this->_weights};
  }

  return iterator<N>{this->_index,         *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                     *this->_pixels_count, this->_coord1,          this->_coord2,
                     this->_weights};
}

template <typename N>
inline auto PixelSelector::cend() const -> iterator<N> {
  if constexpr (std::is_integral_v<N>) {
    if (!!_weights) {
      throw std::logic_error(
          "iterator template parameter should be of floating point type when processing balanced "
          "matrices.");
    }
  }
  return iterator<N>::at_end(this->_index, *this->_pixels_bin1_id, *this->_pixels_bin2_id,
                             *this->_pixels_count, this->_weights);
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all() const {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel<N>> buff{};
  std::transform(begin<N>(), end<N>(), std::back_inserter(buff), [&](const ThinPixel<N> &p) {
    return Pixel<N>{{_index->bins().at(p.bin1_id), _index->bins().at(p.bin2_id)}, p.count};
  });
  return buff;
}

inline const PixelCoordinates &PixelSelector::coord1() const noexcept { return this->_coord1; }

inline const PixelCoordinates &PixelSelector::coord2() const noexcept { return this->_coord2; }

inline const BinTable &PixelSelector::bins() const noexcept { return this->_index->bins(); }

inline std::shared_ptr<const BinTable> PixelSelector::bins_ptr() const noexcept {
  return this->_index->bins_ptr();
}

template <typename N>
inline PixelSelector::iterator<N>::iterator(std::shared_ptr<const Index> index,
                                            const Dataset &pixels_bin1_id,
                                            const Dataset &pixels_bin2_id,
                                            const Dataset &pixels_count,
                                            std::shared_ptr<const balancing::Weights> weights)
    : _bin1_id_it(pixels_bin1_id.begin<BinIDT>(32'000)),
      _bin2_id_it(pixels_bin2_id.begin<BinIDT>(32'000)),
      _count_it(pixels_count.begin<N>(32'000)),
      _index(std::move(index)),
      _weights(std::move(weights)),
      _h5_end_offset(pixels_bin2_id.size()) {
  std::ignore = **this;
}

template <typename N>
inline PixelSelector::iterator<N>::iterator(std::shared_ptr<const Index> index,
                                            const Dataset &pixels_bin1_id,
                                            const Dataset &pixels_bin2_id,
                                            const Dataset &pixels_count, PixelCoordinates coord1,
                                            PixelCoordinates coord2,
                                            std::shared_ptr<const balancing::Weights> weights)
    : _index(std::move(index)),
      _coord1(std::move(coord1)),
      _coord2(std::move(coord2)),
      _weights(std::move(weights)),
      _h5_end_offset(pixels_bin2_id.size()) {
  assert(_coord1);
  assert(_coord2);
  assert(_coord1.bin1.id() <= _coord1.bin2.id());
  assert(_coord2.bin1.id() <= _coord2.bin2.id());

  const auto chunk_size = 32'000;
  // Set iterator to the first row overlapping the query (i.e. the first bin overlapping coord1)
  auto offset = _index->get_offset_by_bin_id(_coord1.bin1.id());
  _bin1_id_it = pixels_bin1_id.make_iterator_at_offset<BinIDT>(offset, chunk_size);
  _bin2_id_it = pixels_bin2_id.make_iterator_at_offset<BinIDT>(offset, chunk_size);
  _count_it = pixels_count.make_iterator_at_offset<N>(offset, chunk_size);

  // Now that last it is set, we can call jump_to_col() to seek to the first pixel actually
  // overlapping the query. Calling jump_to_next_overlap() is required to deal with rows that are
  // not empty, but that have no pixels overlapping the query
  this->jump_to_col(_coord2.bin1.id());
  if (this->discard()) {
    this->jump_to_next_overlap();
  }

  if (this->is_at_end()) {
    *this = at_end(std::move(this->_index), pixels_bin1_id, pixels_bin2_id, pixels_count, _weights);
  } else {
    std::ignore = **this;
  }
}

template <typename N>
inline auto PixelSelector::iterator<N>::at_end(std::shared_ptr<const Index> index,
                                               const Dataset &pixels_bin1_id,
                                               const Dataset &pixels_bin2_id,
                                               const Dataset &pixels_count,
                                               std::shared_ptr<const balancing::Weights> weights)
    -> iterator {
  iterator it{};
  it._index = std::move(index);
  it._bin1_id_it = pixels_bin1_id.end<BinIDT>(0);
  it._bin2_id_it = pixels_bin2_id.end<BinIDT>(0);
  it._count_it = pixels_count.end<N>(0);
  it._weights = std::move(weights);
  it._h5_end_offset = pixels_bin2_id.size();

  return it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator==(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it == other._bin2_id_it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it < other._bin2_id_it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<=(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it <= other._bin2_id_it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator>(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it > other._bin2_id_it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator>=(const iterator &other) const noexcept {
  assert(this->_index == other._index);
  return this->_bin2_id_it >= other._bin2_id_it;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator*() const -> const_reference {
  assert(!this->is_at_end());
  _value = {*this->_bin1_id_it, *this->_bin2_id_it, conditional_static_cast<N>(*this->_count_it)};

  if constexpr (std::is_floating_point_v<N>) {
    if (this->_weights) {
      _value = this->_weights->balance(_value);
    }
  }
  return _value;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator->() const -> const_pointer {
  return &(**this);
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++() -> iterator & {
  assert(!this->is_at_end());
  std::ignore = ++this->_bin1_id_it;
  std::ignore = ++this->_bin2_id_it;
  std::ignore = ++this->_count_it;

  if (this->is_at_end()) {
    this->jump_at_end();
    return *this;
  }

  if (this->discard()) {
    this->jump_to_next_overlap();
  }

  return *this;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++(int) -> iterator {
  if (this->_bin1_id_it.underlying_buff_num_available_fwd() <= 1) {
    this->refresh();
  }
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N>
inline void PixelSelector::iterator<N>::jump_to_row(std::uint64_t bin_id) {
  assert(this->_index);
  assert(bin_id <= this->_index->bins().size());

  if (this->is_at_end()) {
    return;
  }

  const auto row_offset = this->_index->get_offset_by_bin_id(bin_id);
  const auto current_offset = this->h5_offset();

  assert(row_offset >= current_offset);
  const auto offset = row_offset - current_offset;

  this->_bin1_id_it += offset;
  this->_bin2_id_it += offset;
  this->_count_it += offset;
}

template <typename N>
inline void PixelSelector::iterator<N>::jump_to_col(std::uint64_t bin_id) {
  assert(this->_index);
  assert(bin_id <= this->_index->bins().size());

  if (this->is_at_end()) {
    return;
  }

  const auto current_row = *this->_bin1_id_it;
  const auto next_row = current_row + 1;

  const auto current_offset = conditional_static_cast<std::uint64_t>(this->h5_offset());
  const auto current_row_offset = this->_index->get_offset_by_bin_id(current_row);
  const auto next_row_offset = this->_index->get_offset_by_bin_id(next_row);

  if (current_offset == next_row_offset) {
    return;  // Row is empty
  }

  assert(next_row_offset != 0);
  const auto row_start_offset = (std::min)(current_offset, current_row_offset);
  const auto row_end_offset = next_row_offset - 1;

  if (row_start_offset == row_end_offset) {
    return;  // Row is empty
  }

  // Doing binary search on iterators can be expensive.
  // Reading an entire row on bin_ids into a vector,
  // and then doing the binary search on the vector is in practice ~15% faster
  const auto &dset = this->_bin2_id_it.dataset();
  const auto chunk_size = row_end_offset - row_start_offset;
  const auto offset1 = this->_bin2_id_it.h5_offset() - (current_offset - row_start_offset);
  auto first = dset.template make_iterator_at_offset<BinIDT>(offset1, chunk_size);
  auto last = dset.template make_iterator_at_offset<BinIDT>(offset1 + chunk_size, 0);
  this->_bin2_id_it = std::lower_bound(first, last, bin_id);

  const auto offset = this->_bin2_id_it.h5_offset() - current_offset;

  this->_bin1_id_it += offset;
  this->_count_it += offset;

  assert(*this->_bin1_id_it == current_row);
}

template <typename N>
inline void PixelSelector::iterator<N>::jump(std::uint64_t bin1_id, std::uint64_t bin2_id) {
  assert(bin1_id <= bin2_id);

  this->jump_to_row(bin1_id);
  if (bin2_id != bin1_id) {
    this->jump_to_col(bin2_id);
  }
}

template <typename N>
inline void PixelSelector::iterator<N>::jump_to_next_overlap() {
  assert(this->discard());
  assert(this->_coord1);
  assert(this->_coord2);
  do {
    // We're at/past end: return immediately
    if (this->is_at_end()) {
      this->jump_at_end();
      return;
    }

    const auto row = *this->_bin1_id_it;
    const auto col = *this->_bin2_id_it;
    const auto next_row = row + 1;
    const auto next_col = (std::max)(next_row, this->_coord2.bin1.id());

    // We may have some data left to read from the current row
    if (col < this->_coord2.bin1.id()) {
      this->jump_to_col(this->_coord2.bin1.id());
      if (!this->discard()) {
        return;
      }
    }

    // There's no more data to be read, as we're past the last column overlapping the query,
    // and the next row does not overlap the query
    if (this->is_at_end() || next_row > this->_coord1.bin2.id()) {
      // assert(col > this->_coord2.bin2.id());  // This is not always true for trans queries
      this->jump_at_end();
      return;
    }

    this->jump(next_row, next_col);
  } while (this->discard());

  if (this->is_at_end()) {
    this->jump_at_end();
    return;
  }
}

template <typename N>
inline std::size_t PixelSelector::iterator<N>::h5_offset() const noexcept {
  assert(this->_bin1_id_it.h5_offset() == this->_bin2_id_it.h5_offset());
  assert(this->_count_it.h5_offset() == this->_bin2_id_it.h5_offset());

  return this->_bin2_id_it.h5_offset();
}

template <typename N>
inline void PixelSelector::iterator<N>::jump_at_end() {
  if (this->_h5_end_offset != this->_bin2_id_it.h5_offset()) {
    *this = at_end(std::move(this->_index), this->_bin1_id_it.dataset(),
                   this->_bin2_id_it.dataset(), this->_count_it.dataset(), this->_weights);
  }
}

template <typename N>
inline void PixelSelector::iterator<N>::refresh() {
  const auto h5_offset = this->_bin1_id_it.h5_offset();

  const auto &bin1_dset = this->_bin1_id_it.dataset();
  const auto &bin2_dset = this->_bin2_id_it.dataset();
  const auto &count_dset = this->_count_it.dataset();

  const auto chunk_size = 32'000;
  this->_bin1_id_it = bin1_dset.template make_iterator_at_offset<BinIDT>(h5_offset, chunk_size);
  this->_bin2_id_it = bin2_dset.template make_iterator_at_offset<BinIDT>(h5_offset, chunk_size);
  this->_count_it = count_dset.template make_iterator_at_offset<N>(h5_offset, chunk_size);
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::overlaps_coord1() const noexcept {
  return !this->_coord1 || (*this->_bin1_id_it >= this->_coord1.bin1.id() &&
                            *this->_bin1_id_it <= this->_coord1.bin2.id());
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::overlaps_coord2() const noexcept {
  return !this->_coord2 || (*this->_bin2_id_it >= this->_coord2.bin1.id() &&
                            *this->_bin2_id_it <= this->_coord2.bin2.id());
}

template <typename N>
inline bool PixelSelector::iterator<N>::discard() const {
  if (this->is_at_end()) {
    return false;
  }

  return !this->overlaps_coord1() || !this->overlaps_coord2();
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::is_at_end() const noexcept {
  if (this->_h5_end_offset == this->_bin2_id_it.h5_offset()) {
    return true;
  }
  return !this->overlaps_coord1() && !this->overlaps_coord2();
}

}  // namespace hictk::cooler
