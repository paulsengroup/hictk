// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/pixel.hpp"

namespace hictk::cooler {

inline PixelSelector::PixelSelector(const Index &index, const Dataset &pixels_bin1_id,
                                    const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                                    std::shared_ptr<const balancing::Weights> weights,
                                    bool symmetric_upper_) noexcept
    : PixelSelector(index.bins_ptr(), pixels_bin1_id, pixels_bin2_id, pixels_count,
                    std::move(weights), symmetric_upper_) {}

inline PixelSelector::PixelSelector(std::shared_ptr<const BinTable> bins,
                                    const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                    const Dataset &pixels_count,
                                    std::shared_ptr<const balancing::Weights> weights,
                                    bool symmetric_upper_) noexcept
    : _bins(std::move(bins)),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count),
      _weights(std::move(weights)),
      _symmetric_upper(symmetric_upper_) {
  assert(_bins);
  assert(_weights);
}

inline PixelSelector::PixelSelector(std::shared_ptr<const Index> index,
                                    const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                    const Dataset &pixels_count, const PixelCoordinates &coords,
                                    std::shared_ptr<const balancing::Weights> weights,
                                    bool symmetric_upper_)
    : PixelSelector(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count, coords, coords,
                    std::move(weights), symmetric_upper_) {}

inline PixelSelector::PixelSelector(std::shared_ptr<const Index> index,
                                    const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                                    const Dataset &pixels_count, PixelCoordinates coord1,
                                    PixelCoordinates coord2,
                                    std::shared_ptr<const balancing::Weights> weights,
                                    bool symmetric_upper_)
    : _coord1(std::move(coord1)),
      _coord2(std::move(coord2)),
      _index(std::move(index)),
      _bins(_index->bins_ptr()),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count),
      _weights(std::move(weights)),
      _symmetric_upper(symmetric_upper_) {
  assert(_index);
  assert(_weights);

  if (!_symmetric_upper) {
    assert(_coord1.empty() && _coord2.empty());
  }

  const auto query_is_cis = _coord1.bin1.chrom() == _coord2.bin1.chrom();
  if ((!query_is_cis && _coord1.bin1 > _coord2.bin1) ||
      (query_is_cis && _coord1.bin1.start() > _coord2.bin1.start())) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("query {}:{}-{}; {}:{}-{}; overlaps with the lower-triangle of the matrix"),
        _coord1.bin1.chrom().name(), _coord1.bin1.start(), _coord1.bin2.end(),
        _coord2.bin1.chrom().name(), _coord2.bin1.start(), _coord2.bin2.end()));
  }
}

inline bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  // clang-format off
  return begin<int>() == other.begin<int>() &&
         end<int>() == other.end<int>() &&
         _weights == other._weights &&
         _symmetric_upper == other._symmetric_upper;
  // clang-format on
}

inline bool PixelSelector::operator!=(const PixelSelector &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelector::begin() const -> iterator<N> {
  return cbegin<N>();
}

template <typename N>
inline auto PixelSelector::end() const -> iterator<N> {
  return cend<N>();
}

template <typename N>
inline auto PixelSelector::cbegin() const -> iterator<N> {
  if constexpr (std::is_integral_v<N>) {
    assert(_weights);
    if (!_weights->is_vector_of_ones()) {
      throw std::logic_error(
          "iterator template parameter should be of floating point type when processing balanced "
          "matrices.");
    }
  }

  const auto fixed_bins = _bins->type() == BinTable::Type::fixed;

  if (!_coord1) {
    assert(!_coord2);
    return iterator<N>{*_pixels_bin1_id, *_pixels_bin2_id, *_pixels_count,
                       _weights,         fixed_bins,       _index};
  }

  return iterator<N>{_index,  *_pixels_bin1_id, *_pixels_bin2_id, *_pixels_count,
                     _coord1, _coord2,          _weights,         fixed_bins};
}

template <typename N>
inline auto PixelSelector::cend() const -> iterator<N> {
  if constexpr (std::is_integral_v<N>) {
    assert(_weights);
    if (!_weights->is_vector_of_ones()) {
      throw std::logic_error(
          "iterator template parameter should be of floating point type when processing balanced "
          "matrices.");
    }
  }

  const auto fixed_bins = _bins->type() == BinTable::Type::fixed;

  return iterator<N>::at_end(_index, *_pixels_bin1_id, *_pixels_bin2_id, *_pixels_count, _weights,
                             fixed_bins);
}

inline bool PixelSelector::empty() const { return begin<double>() == end<double>(); }

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all() const {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel<N>> buff{};
  std::transform(begin<N>(), end<N>(), std::back_inserter(buff), [&](const ThinPixel<N> &p) {
    return Pixel<N>{{_bins->at(p.bin1_id), _bins->at(p.bin2_id)}, p.count};
  });
  return buff;
}

inline const PixelCoordinates &PixelSelector::coord1() const noexcept { return _coord1; }

inline const PixelCoordinates &PixelSelector::coord2() const noexcept { return _coord2; }

inline const BinTable &PixelSelector::bins() const noexcept { return *bins_ptr(); }

inline std::shared_ptr<const BinTable> PixelSelector::bins_ptr() const noexcept { return _bins; }

inline PixelSelector PixelSelector::fetch(PixelCoordinates coord1, PixelCoordinates coord2) const {
  return {_index,         *_pixels_bin1_id,  *_pixels_bin2_id,
          *_pixels_count, std::move(coord1), std::move(coord2),
          _weights,       _symmetric_upper};
}

inline const balancing::Weights &PixelSelector::weights() const noexcept {
  assert(_weights);
  return *_weights;
}

inline bool PixelSelector::is_symmetric_upper() const noexcept { return _symmetric_upper; }

template <typename N>
inline PixelSelector::iterator<N>::iterator(
    // NOLINTBEGIN(*-unnecessary-value-param)
    const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id, const Dataset &pixels_count,
    std::shared_ptr<const balancing::Weights> weights, bool fixed_bin_size,
    std::shared_ptr<const Index> index)
    // NOLINTEND(*-unnecessary-value-param)
    : _bin1_id_it(pixels_bin1_id.begin<BinIDT>()),
      _bin2_id_it(pixels_bin2_id.begin<BinIDT>()),
      _count_it(pixels_count.begin<N>()),
      _index(std::move(index)),
      _weights(std::move(weights)),
      _h5_end_offset(pixels_bin2_id.size()),
      _fixed_bin_size(fixed_bin_size) {
  assert(_weights);
}

template <typename N>
inline PixelSelector::iterator<N>::iterator(
    // NOLINTBEGIN(*-unnecessary-value-param)
    std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
    const Dataset &pixels_bin2_id, const Dataset &pixels_count, PixelCoordinates coord1,
    PixelCoordinates coord2, std::shared_ptr<const balancing::Weights> weights, bool fixed_bin_size)
    // NOLINTEND(*-unnecessary-value-param)
    : _index(std::move(index)),
      _coord1(std::move(coord1)),
      _coord2(std::move(coord2)),
      _weights(std::move(weights)),
      _h5_end_offset(pixels_bin2_id.size()),
      _fixed_bin_size(fixed_bin_size) {
  assert(_coord1);
  assert(_coord2);
  assert(_coord1.bin1.id() <= _coord1.bin2.id());
  assert(_coord2.bin1.id() <= _coord2.bin2.id());
  assert(_weights);

  if (_index->empty(_coord1.bin1.chrom().id())) {
    *this = at_end(std::move(_index), pixels_bin1_id, pixels_bin2_id, pixels_count,
                   std::move(_weights), _fixed_bin_size);
    return;
  }

  // Set iterator to the first row overlapping the query (i.e. the first bin overlapping coord1)
  const auto offset = _index->get_offset_by_bin_id(_coord1.bin1.id());
  _bin1_id_it = pixels_bin1_id.make_iterator_at_offset<BinIDT>(offset);
  _bin2_id_it = pixels_bin2_id.make_iterator_at_offset<BinIDT>(offset);
  _count_it = pixels_count.make_iterator_at_offset<N>(offset);

  // Now that last it is set, we can call jump_to_col() to seek to the first pixel actually
  // overlapping the query. Calling jump_to_next_overlap() is required to deal with rows that are
  // not empty, but that have no pixels overlapping the query
  jump_to_col(_coord2.bin1.id());
  if (discard()) {
    jump_to_next_overlap();
  }

  if (is_at_end()) {
    *this = at_end(std::move(_index), pixels_bin1_id, pixels_bin2_id, pixels_count, _weights,
                   _fixed_bin_size);
  }
}

template <typename N>
inline auto PixelSelector::iterator<N>::at_end(std::shared_ptr<const Index> index,
                                               const Dataset &pixels_bin1_id,
                                               const Dataset &pixels_bin2_id,
                                               const Dataset &pixels_count,
                                               std::shared_ptr<const balancing::Weights> weights,
                                               bool fixed_bin_size) -> iterator {
  iterator it{};
  it._index = std::move(index);
  it._bin1_id_it = pixels_bin1_id.end<BinIDT>(0);
  it._bin2_id_it = pixels_bin2_id.end<BinIDT>(0);
  it._count_it = pixels_count.end<N>(0);
  it._weights = std::move(weights);
  it._h5_end_offset = pixels_bin2_id.size();
  it._fixed_bin_size = fixed_bin_size;

  return it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator==(const iterator &other) const noexcept {
  assert(_index == other._index);
  return _bin2_id_it == other._bin2_id_it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<(const iterator &other) const noexcept {
  assert(_index == other._index);
  return _bin2_id_it < other._bin2_id_it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<=(const iterator &other) const noexcept {
  assert(_index == other._index);
  return _bin2_id_it <= other._bin2_id_it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator>(const iterator &other) const noexcept {
  assert(_index == other._index);
  return _bin2_id_it > other._bin2_id_it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator>=(const iterator &other) const noexcept {
  assert(_index == other._index);
  return _bin2_id_it >= other._bin2_id_it;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator*() const -> const_reference {
  assert(!is_at_end());
  _value = {*_bin1_id_it, *_bin2_id_it, conditional_static_cast<N>(*_count_it)};

  if constexpr (std::is_floating_point_v<N>) {
    assert(_weights);
    _value.count = _weights->balance<N>(_value.bin1_id, _value.bin2_id, _value.count);
  } else {
    assert(_weights);
    assert(_weights->is_vector_of_ones());
  }
  return _value;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator->() const -> const_pointer {
  return &(**this);
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++() -> iterator & {
  assert(!is_at_end());
  std::ignore = ++_bin1_id_it;
  std::ignore = ++_bin2_id_it;
  std::ignore = ++_count_it;

  if (is_at_end()) {
    jump_at_end();
    return *this;
  }

  if (discard()) {
    jump_to_next_overlap();
  }

  return *this;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++(int) -> iterator {
  if (_bin1_id_it.h5_offset() + 1 >= _bin1_id_it.buffer().end()) {
    refresh();
  }
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N>
inline void PixelSelector::iterator<N>::jump_to_row(std::uint64_t bin_id) {
  assert(_index);
  assert(bin_id <= _index->bins().size());

  if (is_at_end()) {
    return;
  }

  const auto row_offset = _index->get_offset_by_bin_id(bin_id);
  const auto current_offset = h5_offset();

  assert(row_offset >= current_offset);
  const auto offset = static_cast<std::ptrdiff_t>(row_offset - current_offset);

  _bin1_id_it += offset;
  _bin2_id_it += offset;
  _count_it += offset;
}

template <typename N>
inline void PixelSelector::iterator<N>::jump_to_col(std::uint64_t bin_id) {
  assert(_index);
  assert(bin_id <= _index->bins().size());

  if (is_at_end()) {
    return;
  }

  const auto current_row = *_bin1_id_it;
  const auto next_row = current_row + 1;

  const auto current_offset = conditional_static_cast<std::uint64_t>(h5_offset());
  const auto current_row_offset = _index->get_offset_by_bin_id(current_row);
  const auto next_row_offset = _index->get_offset_by_bin_id(next_row);

  if (current_offset == next_row_offset) {
    return;  // Row is empty
  }

  assert(next_row_offset != 0);
  const auto row_start_offset = (std::min)(current_offset, current_row_offset);
  const auto row_end_offset = next_row_offset - 1;

  if (row_start_offset == row_end_offset) {
    return;  // Row is empty
  }

  assert(row_end_offset != 0);
  const auto row_size = row_end_offset - row_start_offset;
  const auto offset = row_start_offset;

  const auto optimistic_offset =
      std::clamp(offset + _row_head_h5_offset, row_start_offset, row_end_offset - 1);
  auto first = _bin2_id_it.seek(optimistic_offset);
  if (HICTK_LIKELY(_coord1 == _coord2 && *first >= bin_id)) {
    // this is mostly an optimization for symmetric cis queries
    assert(optimistic_offset == offset);
    _bin2_id_it = std::move(first);
  } else if (HICTK_LIKELY(*first <= bin_id)) {
    // optimistic first is at/upstream of the bin_id
    auto &last = _bin2_id_it.seek(offset + row_size);
    _bin2_id_it = Dataset::lower_bound(std::move(first), std::move(last), bin_id, true);
  } else {
    // optimistic first is downstream of the bin_id: need to search through the entire row
    if (optimistic_offset != offset) {
      first.seek(offset);
    }
    auto &last = _bin2_id_it.seek(offset + row_size);
    _bin2_id_it = Dataset::lower_bound(std::move(first), std::move(last), bin_id, true);
  }

  _bin1_id_it.seek(_bin2_id_it.h5_offset());
  _count_it.seek(_bin2_id_it.h5_offset());

  if (_coord1 != _coord2) {
    _row_head_h5_offset = h5_offset() - row_start_offset;

    // try to reduce the offset by 5% of the row length to increase the likelihood that
    // row_start_offset + _row_head_h5_offset is upstream of the given bin_id
    const auto left_shift_offset = (5 * row_size) / 100;
    _row_head_h5_offset -= std::min(_row_head_h5_offset, left_shift_offset);
  }

  assert(*_bin1_id_it == current_row);
}

template <typename N>
inline void PixelSelector::iterator<N>::jump(std::uint64_t bin1_id, std::uint64_t bin2_id) {
  assert(bin1_id <= bin2_id);

  jump_to_row(bin1_id);
  if (bin2_id != bin1_id) {
    jump_to_col(bin2_id);
  }
}

template <typename N>
inline bool PixelSelector::iterator<N>::is_indexed() const noexcept {
  if (!_index) {
    return false;
  }

  if (!_coord1 && !_coord2) {
    return std::all_of(_index->chromosomes().begin(), _index->chromosomes().end(),
                       [&](const Chromosome &chrom) { return _index->contains(chrom.id()); });
  }

  bool loaded = true;
  if (_coord1) {
    loaded &= _index->contains(_coord1.bin1.chrom().id());
  }

  if (_coord2) {
    loaded &= _index->contains(_coord2.bin1.chrom().id());
  }

  return loaded;
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::is_fixed_bin_size() const noexcept {
  return _fixed_bin_size;
}

template <typename N>
inline void PixelSelector::iterator<N>::jump_to_next_overlap() {
  if (!_coord1) {
    assert(!_coord2);
    return;
  }

  do {
    // We're at/past end: return immediately
    if (is_at_end()) {
      jump_at_end();
      return;
    }

    const auto row = *_bin1_id_it;
    const auto col = *_bin2_id_it;
    const auto next_row = row + 1;
    const auto next_col = (std::max)(next_row, _coord2.bin1.id());

    // We may have some data left to read from the current row
    if (col < _coord2.bin1.id()) {
      jump_to_col(_coord2.bin1.id());
      if (!discard()) {
        return;
      }
    }

    // There's no more data to be read, as we're past the last column overlapping the query,
    // and the next row does not overlap the query
    if (is_at_end() || next_row > _coord1.bin2.id()) {
      // assert(col > _coord2.bin2.id());  // This is not always true for trans queries
      jump_at_end();
      return;
    }

    jump(next_row, next_col);
  } while (discard());

  if (is_at_end()) {
    jump_at_end();
  }
}

template <typename N>
inline std::size_t PixelSelector::iterator<N>::h5_offset() const noexcept {
  assert(_bin1_id_it.h5_offset() == _bin2_id_it.h5_offset());
  assert(_count_it.h5_offset() == _bin2_id_it.h5_offset());

  return _bin2_id_it.h5_offset();
}

template <typename N>
inline void PixelSelector::iterator<N>::jump_at_end() {
  if (_h5_end_offset != _bin2_id_it.h5_offset()) {
    *this = at_end(std::move(_index), _bin1_id_it.dataset(), _bin2_id_it.dataset(),
                   _count_it.dataset(), _weights, _fixed_bin_size);
  }
}

template <typename N>
inline void PixelSelector::iterator<N>::refresh() {
  const auto h5_offset = _bin1_id_it.h5_offset();

  const auto &bin1_dset = _bin1_id_it.dataset();
  const auto &bin2_dset = _bin2_id_it.dataset();
  const auto &count_dset = _count_it.dataset();

  _bin1_id_it = bin1_dset.template make_iterator_at_offset<BinIDT>(h5_offset,
                                                                   _bin1_id_it.buffer().capacity());
  _bin2_id_it = bin2_dset.template make_iterator_at_offset<BinIDT>(h5_offset,
                                                                   _bin2_id_it.buffer().capacity());
  _count_it =
      count_dset.template make_iterator_at_offset<N>(h5_offset, _count_it.buffer().capacity());
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::overlaps_coord1() const {
  return !_coord1 || (*_bin1_id_it >= _coord1.bin1.id() && *_bin1_id_it <= _coord1.bin2.id());
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::overlaps_coord2() const {
  return !_coord2 || (*_bin2_id_it >= _coord2.bin1.id() && *_bin2_id_it <= _coord2.bin2.id());
}

template <typename N>
inline bool PixelSelector::iterator<N>::discard() const {
  if (is_at_end()) {
    return false;
  }

  return !overlaps_coord1() || !overlaps_coord2();
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::is_at_end() const {
  if (_h5_end_offset == _bin2_id_it.h5_offset()) {
    return true;
  }
  return !overlaps_coord1() && !overlaps_coord2();
}

}  // namespace hictk::cooler
