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
      _bins(_index->bins_ptr()),
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
    : _bins(index->bins_ptr()),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count),
      _weights(std::move(weights)) {}

inline bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  // clang-format off
  return begin<int>() == other.begin<int>() &&
         end<int>() == other.end<int>() &&
         _weights == other._weights;
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
    if (!!_weights) {
      throw std::logic_error(
          "iterator template parameter should be of floating point type when processing balanced "
          "matrices.");
    }
  }

  if (!_coord1) {
    assert(!_coord2);
    return iterator<N>{*_pixels_bin1_id, *_pixels_bin2_id, *_pixels_count, _weights};
  }

  return iterator<N>{_index,  *_pixels_bin1_id, *_pixels_bin2_id, *_pixels_count,
                     _coord1, _coord2,          _weights};
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
  return iterator<N>::at_end(_index, *_pixels_bin1_id, *_pixels_bin2_id, *_pixels_count, _weights);
}

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

#ifdef HICTK_WITH_EIGEN
template <typename N>
inline Eigen::SparseMatrix<N> PixelSelector::read_sparse() const {
  const auto bin_size = _bins->bin_size();
  const auto span1 = coord1().bin2.end() - coord1().bin1.start();
  const auto span2 = coord2().bin2.end() - coord2().bin1.start();
  const auto num_rows = static_cast<std::int64_t>((span1 + bin_size - 1) / bin_size);
  const auto num_cols = static_cast<std::int64_t>((span2 + bin_size - 1) / bin_size);

  const auto offset1 = coord1().bin1.id();
  const auto offset2 = coord2().bin1.id();

  Eigen::SparseMatrix<N> matrix(num_rows, num_cols);
  std::for_each(begin<N>(), end<N>(), [&](const ThinPixel<N> &p) {
    matrix.insert(static_cast<std::int64_t>(p.bin1_id - offset1),
                  static_cast<std::int64_t>(p.bin2_id - offset2)) = p.count;
  });
  matrix.makeCompressed();
  return matrix;
}

template <typename N>
[[nodiscard]] Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic> PixelSelector::read_dense() const {
  const auto bin_size = _bins->bin_size();
  const auto span1 = coord1().bin2.end() - coord1().bin1.start();
  const auto span2 = coord2().bin2.end() - coord2().bin1.start();
  const auto num_rows = static_cast<std::int64_t>((span1 + bin_size - 1) / bin_size);
  const auto num_cols = static_cast<std::int64_t>((span2 + bin_size - 1) / bin_size);

  const auto offset1 = coord1().bin1.id();
  const auto offset2 = coord2().bin1.id();

  using MatrixT = Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  MatrixT matrix = MatrixT::Zero(num_rows, num_cols);
  std::for_each(begin<N>(), end<N>(), [&](const ThinPixel<N> &p) {
    const auto i1 = static_cast<std::int64_t>(p.bin1_id - offset1);
    const auto i2 = static_cast<std::int64_t>(p.bin2_id - offset2);
    matrix(i1, i2) = p.count;

    //  Mirror matrix below diagonal
    if (i2 - i1 < num_rows && i1 < num_cols && i2 < num_rows) {
      matrix(i2, i1) = p.count;
    } else if (i2 - i1 > num_cols && i1 < num_cols && i2 < num_rows) {
      const auto i3 = static_cast<std::int64_t>(p.bin2_id - offset1);
      const auto i4 = static_cast<std::int64_t>(p.bin1_id - offset2);
      matrix(i3, i4) = p.count;
    }
  });
  return matrix;
}
#endif

inline const PixelCoordinates &PixelSelector::coord1() const noexcept { return _coord1; }

inline const PixelCoordinates &PixelSelector::coord2() const noexcept { return _coord2; }

inline const BinTable &PixelSelector::bins() const noexcept { return *bins_ptr(); }

inline std::shared_ptr<const BinTable> PixelSelector::bins_ptr() const noexcept { return _bins; }

template <typename N>
inline PixelSelector::iterator<N>::iterator(const Dataset &pixels_bin1_id,
                                            const Dataset &pixels_bin2_id,
                                            const Dataset &pixels_count,
                                            std::shared_ptr<const balancing::Weights> weights)
    : _bin1_id_it(pixels_bin1_id.begin<BinIDT>()),
      _bin2_id_it(pixels_bin2_id.begin<BinIDT>()),
      _count_it(pixels_count.begin<N>()),
      _weights(std::move(weights)),
      _h5_end_offset(pixels_bin2_id.size()) {}

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

  // Set iterator to the first row overlapping the query (i.e. the first bin overlapping coord1)
  auto offset = _index->get_offset_by_bin_id(_coord1.bin1.id());
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
    *this = at_end(std::move(_index), pixels_bin1_id, pixels_bin2_id, pixels_count, _weights);
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
    if (_weights) {
      _value = _weights->balance(_value);
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
  if (_bin1_id_it.underlying_buff_num_available_fwd() <= 1) {
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
  const auto offset = row_offset - current_offset;

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

  const auto chunk_size = row_end_offset - row_start_offset;
  const auto offset = _bin2_id_it.h5_offset() - (current_offset - row_start_offset);
  auto first = _bin2_id_it.seek(offset);
  auto last = _bin2_id_it.seek(offset + chunk_size);
  _bin2_id_it = std::lower_bound(first, last, bin_id);

  _bin1_id_it.seek(_bin2_id_it.h5_offset());
  _count_it.seek(_bin2_id_it.h5_offset());

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
inline void PixelSelector::iterator<N>::jump_to_next_overlap() {
  assert(discard());
  assert(_coord1);
  assert(_coord2);
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
    return;
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
                   _count_it.dataset(), _weights);
  }
}

template <typename N>
inline void PixelSelector::iterator<N>::refresh() {
  const auto h5_offset = _bin1_id_it.h5_offset();

  const auto &bin1_dset = _bin1_id_it.dataset();
  const auto &bin2_dset = _bin2_id_it.dataset();
  const auto &count_dset = _count_it.dataset();

  _bin1_id_it = bin1_dset.template make_iterator_at_offset<BinIDT>(h5_offset);
  _bin2_id_it = bin2_dset.template make_iterator_at_offset<BinIDT>(h5_offset);
  _count_it = count_dset.template make_iterator_at_offset<N>(h5_offset);
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::overlaps_coord1() const noexcept {
  return !_coord1 || (*_bin1_id_it >= _coord1.bin1.id() && *_bin1_id_it <= _coord1.bin2.id());
}

template <typename N>
constexpr bool PixelSelector::iterator<N>::overlaps_coord2() const noexcept {
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
constexpr bool PixelSelector::iterator<N>::is_at_end() const noexcept {
  if (_h5_end_offset == _bin2_id_it.h5_offset()) {
    return true;
  }
  return !overlaps_coord1() && !overlaps_coord2();
}

}  // namespace hictk::cooler
