// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <random>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/hic/block_cache.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic {

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockCache> cache_,
                                    std::shared_ptr<const BinTable> bins_,
                                    PixelCoordinates coords) noexcept
    : PixelSelector(std::move(hfs_), std::move(footer_), std::move(cache_), std::move(bins_),
                    coords, std::move(coords)) {}

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockCache> cache_,
                                    std::shared_ptr<const BinTable> bins_, PixelCoordinates coord1_,
                                    PixelCoordinates coord2_) noexcept
    : _reader(std::move(hfs_), footer_->index(), std::move(bins_), std::move(cache_)),
      _footer(std::move(footer_)),
      _coord1(std::move(coord1_)),
      _coord2(std::move(coord2_)) {}

inline bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  return _reader.index().chrom1() == _reader.index().chrom2() && _coord1 == other._coord1 &&
         _coord2 == other._coord2;
}

inline bool PixelSelector::operator!=(const PixelSelector &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelector::cbegin() const -> iterator<N> {
  return iterator<N>(*this);
}

template <typename N>
inline auto PixelSelector::cend() const -> iterator<N> {
  return iterator<N>::at_end(*this);
}

template <typename N>
inline auto PixelSelector::begin() const -> iterator<N> {
  return this->cbegin<N>();
}

template <typename N>
inline auto PixelSelector::end() const -> iterator<N> {
  return this->cend<N>();
}

inline SerializedPixel PixelSelector::transform_pixel(SerializedPixel pixel) const {
  const auto &c1Norm = _footer->c1Norm();
  const auto &c2Norm = _footer->c2Norm();
  const auto &expected = _footer->expectedValues();

  const auto bin1 = static_cast<std::size_t>(pixel.bin1_id);
  const auto bin2 = static_cast<std::size_t>(pixel.bin2_id);

  assert(is_inter() || bin1 <= bin2);

  const auto skipNormalization =
      normalization() == NormalizationMethod::NONE || matrix_type() == MatrixType::expected;

  if (!skipNormalization) {
    assert(bin1 < c1Norm.size());
    assert(bin2 < c2Norm.size());
    pixel.count /= static_cast<float>(c1Norm[bin1] * c2Norm[bin2]);
  }

  if (matrix_type() == MatrixType::observed) {
    return pixel;
  }

  const auto expectedCount = [&]() {
    if (is_inter()) {
      return float(_reader.avg());
    }

    const auto i = (bin2 - bin1);
    assert(i < expected.size());
    return float(expected[i]);
  }();

  if (matrix_type() == MatrixType::expected) {
    pixel.count = expectedCount;
    return pixel;
  }

  assert(matrix_type() == MatrixType::oe);
  pixel.count /= expectedCount;

  return pixel;
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all() const {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel<N>> buff{};
  std::copy(begin<N>(), end<N>(), std::back_inserter(buff));
  return buff;
}

inline const PixelCoordinates &PixelSelector::coord1() const noexcept { return _coord1; }
inline const PixelCoordinates &PixelSelector::coord2() const noexcept { return _coord2; }
inline MatrixType PixelSelector::matrix_type() const noexcept { return metadata().matrix_type; }
inline NormalizationMethod PixelSelector::normalization() const noexcept {
  return metadata().normalization;
}
inline MatrixUnit PixelSelector::unit() const noexcept { return _reader.index().unit(); }
inline std::uint32_t PixelSelector::resolution() const noexcept {
  return _reader.index().resolution();
}

inline const Chromosome &PixelSelector::chrom1() const noexcept { return _coord1.bin1.chrom(); }
inline const Chromosome &PixelSelector::chrom2() const noexcept { return _coord2.bin1.chrom(); }

inline const std::vector<double> &PixelSelector::chrom1_norm() const noexcept {
  return _footer->c1Norm();
}
inline const std::vector<double> &PixelSelector::chrom2_norm() const noexcept {
  return _footer->c2Norm();
}

inline const BinTable &PixelSelector::bins() const noexcept { return _reader.bins(); }

inline const internal::HiCFooterMetadata &PixelSelector::metadata() const noexcept {
  assert(!!this->_footer);
  return this->_footer->metadata();
}

inline bool PixelSelector::is_intra() const noexcept { return chrom1() == chrom2(); }

inline bool PixelSelector::is_inter() const noexcept { return !is_intra(); }

template <typename N>
inline N PixelSelector::sum() const noexcept {
  return _reader.sum();
}
inline double PixelSelector::avg() const noexcept { return _reader.avg(); }

inline std::size_t PixelSelector::estimate_optimal_cache_size(std::size_t num_samples) const {
  if (_reader.index().empty()) {
    return 0;  // should we throw instead?
  }

  std::seed_seq sseq({_reader.index().size()});
  std::mt19937_64 rand_eng(sseq);

  // Try to guess the average block size
  std::size_t max_block_size = 0;

  std::vector<std::size_t> block_sizes{};
  std::size_t samples = 0;
  // index is backed by a hashmap, so iteration should be somewhat random
  for (const auto &idx : _reader.index()) {
    auto blk = _reader.read(chrom1(), chrom2(), idx);
    if (blk) {
      samples++;
      max_block_size = (std::max)(blk->size(), max_block_size);
      _reader.evict(*blk);
    }
    if (samples == num_samples) {
      break;
    }
  }

  // Try to guess how many blocks overlap a single row of pixels
  std::size_t max_blocks_per_row = 0;
  const auto bin_size = bins().bin_size();

  const std::size_t first_bin_id = 0;
  const std::size_t last_bin_id =
      bins().at(coord1().bin1.chrom(), coord1().bin1.chrom().size()).rel_id() - 1;
  samples = (std::min)(num_samples, bins().subset(coord1().bin1.chrom()).size());
  for (std::size_t i = 0; i < samples; ++i) {
    const auto bin_id =
        std::uniform_int_distribution<std::size_t>{first_bin_id, last_bin_id}(rand_eng);

    const auto pos1 = static_cast<std::uint32_t>(bin_id * bin_size);
    const auto pos2 = (std::min)(pos1 + bin_size, coord1().bin1.chrom().size());

    const auto coord1_ = PixelCoordinates(bins().at(coord1().bin1.chrom(), pos1),
                                          bins().at(coord1().bin1.chrom(), pos2));

    std::vector<internal::BlockIndex> buffer{};
    _reader.index().find_overlaps(coord1_, coord2(), buffer);
    max_blocks_per_row = (std::max)(max_blocks_per_row, buffer.size());
  }

  return (std::max)(10'000'000UL, max_blocks_per_row * max_block_size);
}

inline void PixelSelector::evict_blocks_from_cache() const {
  std::vector<internal::BlockIndex> buff{};
  _reader.index().find_overlaps(coord1(), coord2(), buff);
  for (const auto &idx : buff) {
    _reader.evict(chrom1(), chrom2(), idx);
  }
}

template <typename N>
inline PixelSelector::iterator<N>::iterator(const PixelSelector &sel)
    : _sel(&sel),
      _bin1_id(coord1().bin1.rel_id()),
      _block_idx_buffer(std::make_shared<BlockIdxBufferT>()),
      _buffer(std::make_shared<BufferT>()) {
  if (_sel->_reader.index().empty()) {
    *this = at_end(sel);
    return;
  }

  while (!!_buffer && _buffer->empty()) {
    read_next_chunk();
  }
}

template <typename N>
inline auto PixelSelector::iterator<N>::at_end(const PixelSelector &sel) -> iterator<N> {
  iterator it{};

  it._sel = &sel;
  it._buffer = nullptr;  // end of queue

  return it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator==(const iterator &other) const noexcept {
  if (_sel != other._sel) {
    return false;
  }

  const auto bin11 = bin1_id();
  const auto bin12 = bin2_id();
  const auto bin21 = other.bin1_id();
  const auto bin22 = other.bin2_id();

  return bin11 == bin21 && bin12 == bin22;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<(const iterator &other) const noexcept {
  const auto bin11 = bin1_id();
  const auto bin21 = other.bin1_id();

  if (bin11 != bin21) {
    return bin11 < bin21;
  }

  const auto bin12 = bin2_id();
  const auto bin22 = other.bin2_id();

  return bin12 < bin22;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator>(const iterator &other) const noexcept {
  const auto bin11 = bin1_id();
  const auto bin21 = other.bin1_id();

  if (bin11 != bin21) {
    return bin11 > bin21;
  }

  const auto bin12 = bin2_id();
  const auto bin22 = other.bin2_id();

  return bin12 > bin22;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator*() const -> const_reference {
  assert(!!_buffer);
  assert(_buffer_i < _buffer->size());
  return (*_buffer)[_buffer_i];
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator->() const -> const_pointer {
  assert(!!_buffer);
  assert(_buffer_i < _buffer->size());
  return &(*_buffer)[_buffer_i];
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++() -> iterator & {
  assert(!!_buffer);

  ++_pixels_processed;
  ++_buffer_i;
  while (!is_at_end() && _buffer_i >= size()) {
    read_next_chunk();
  }

  return *this;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::is_at_end() const noexcept {
  return _buffer == nullptr;
}

template <typename N>
inline const BinTable &PixelSelector::iterator<N>::bins() const noexcept {
  assert(!!_sel);
  return _sel->bins();
}
template <typename N>
inline const PixelCoordinates &PixelSelector::iterator<N>::coord1() const noexcept {
  assert(!!_sel);
  return _sel->coord1();
}
template <typename N>
inline const PixelCoordinates &PixelSelector::iterator<N>::coord2() const noexcept {
  assert(!!_sel);
  return _sel->coord2();
}

template <typename N>
inline std::size_t PixelSelector::iterator<N>::size() const noexcept {
  return !_buffer ? 0 : _buffer->size();
}

template <typename N>
inline std::size_t PixelSelector::iterator<N>::bin1_id() const noexcept {
  return !is_at_end() ? (*this)->coords.bin1.id() : std::numeric_limits<std::size_t>::max();
}
template <typename N>
inline std::size_t PixelSelector::iterator<N>::bin2_id() const noexcept {
  return !is_at_end() ? (*this)->coords.bin2.id() : std::numeric_limits<std::size_t>::max();
}

template <typename N>
inline std::size_t PixelSelector::iterator<N>::compute_chunk_size(double fraction) const noexcept {
  const auto bin_size = bins().bin_size();
  const auto num_bins = (coord1().bin1.chrom().size() + bin_size - 1) / bin_size;
  const auto max_num_bins =
      (std::max)(1U, static_cast<std::uint32_t>(fraction * static_cast<double>(num_bins)));

  const auto end_pos = coord1().bin2.start();
  const auto pos1 = (std::min)(end_pos, static_cast<std::uint32_t>(_bin1_id) * bins().bin_size());
  const auto pos2 = (std::min)(end_pos, pos1 + (max_num_bins * bin_size));

  return (pos2 - pos1 + bin_size - 1) / bin_size;
}

template <typename N>
inline const std::vector<internal::BlockIndex> &
PixelSelector::iterator<N>::find_blocks_overlapping_next_chunk(std::size_t num_bins) {
  const auto bin_size = bins().bin_size();

  const auto end_pos = coord1().bin2.start();
  const auto pos1 = (std::min)(end_pos, static_cast<std::uint32_t>(_bin1_id) * bins().bin_size());
  const auto pos2 = (std::min)(end_pos, pos1 + static_cast<std::uint32_t>((num_bins * bin_size)));

  const auto coord1_ = PixelCoordinates(bins().at(coord1().bin1.chrom(), pos1),
                                        bins().at(coord1().bin1.chrom(), pos2));

  if (_block_idx_buffer.use_count() != 1) {
    _block_idx_buffer = std::make_shared<BlockIdxBufferT>();
  }

  _sel->_reader.index().find_overlaps(coord1_, coord2(), *_block_idx_buffer);
  return *_block_idx_buffer;
}

template <typename N>
inline void PixelSelector::iterator<N>::read_next_chunk() {
  assert(!!_sel);

  if (_bin1_id > coord1().bin2.rel_id()) {
    *this = at_end(*_sel);
    return;
  }

  if (_buffer.use_count() != 1) {
    _buffer = std::make_shared<BufferT>(_buffer->capacity());
  }
  _buffer->clear();
  _buffer_i = 0;

  const auto chunk_size = compute_chunk_size();
  const auto bin_size = bins().bin_size();
  const auto bin1_id_last = _bin1_id + chunk_size;

  const auto blocks = find_blocks_overlapping_next_chunk(chunk_size);
  if (blocks.empty()) {
    _bin1_id = bin1_id_last + 1;
    return;
  }

  for (const auto &block_idx : blocks) {
    const auto block = _sel->_reader.read(coord1().bin1.chrom(), coord2().bin1.chrom(), block_idx);
    auto first = std::lower_bound(block->begin(), block->end(), _bin1_id,
                                  [](const SerializedPixel &pixel, std::size_t bin_id) {
                                    return pixel.bin1_id < static_cast<std::int64_t>(bin_id);
                                  });
    auto last = std::lower_bound(first, block->end(), bin1_id_last + 1,
                                 [](const SerializedPixel &pixel, std::size_t bin_id) {
                                   return pixel.bin1_id < static_cast<std::int64_t>(bin_id);
                                 });

    const auto buffer_size = _buffer->size();
    while (first != last) {
      const auto p = _sel->transform_pixel(*first++);
      if (p.bin2_id < coord2().bin1.rel_id() || p.bin2_id > coord2().bin2.rel_id()) {
        continue;
      }

      const auto pos1 = static_cast<std::uint32_t>(p.bin1_id) * bin_size;
      const auto pos2 = static_cast<std::uint32_t>(p.bin2_id) * bin_size;
      auto coords = PixelCoordinates{bins().at(coord1().bin1.chrom(), pos1),
                                     bins().at(coord2().bin1.chrom(), pos2)};
      if constexpr (std::is_integral_v<N>) {
        _buffer->emplace_back(Pixel<N>{coords, conditional_static_cast<N>(std::round(p.count))});
      } else {
        _buffer->emplace_back(Pixel<N>{coords, conditional_static_cast<N>(p.count)});
      }
    }

    auto sorted_first = _buffer->begin();
    auto sorted_last = sorted_first + static_cast<std::ptrdiff_t>(buffer_size);
    std::inplace_merge(sorted_first, sorted_last, _buffer->end());
  }
  assert(std::is_sorted(_buffer->begin(), _buffer->end()));
  _bin1_id = bin1_id_last + 1;
}

inline PixelSelectorAll::PixelSelectorAll(std::vector<PixelSelector> selectors_) noexcept
    : _selectors(std::move(selectors_)) {}

template <typename N>
inline auto PixelSelectorAll::begin() const -> iterator<N> {
  return cbegin<N>();
}
template <typename N>
inline auto PixelSelectorAll::cbegin() const -> iterator<N> {
  return iterator<N>(*this);
}

template <typename N>
inline auto PixelSelectorAll::end() const -> iterator<N> {
  return cend<N>();
}
template <typename N>
inline auto PixelSelectorAll::cend() const -> iterator<N> {
  return iterator<N>{};
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelectorAll::read_all() const {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel<N>> buff{};
  std::copy(begin<N>(), end<N>(), std::back_inserter(buff));
  return buff;
}

inline MatrixType PixelSelectorAll::matrix_type() const noexcept {
  return _selectors.front().matrix_type();
}
inline NormalizationMethod PixelSelectorAll::normalization() const noexcept {
  return _selectors.front().normalization();
}
inline MatrixUnit PixelSelectorAll::unit() const noexcept { return _selectors.front().unit(); }
inline std::uint32_t PixelSelectorAll::resolution() const noexcept {
  return _selectors.front().resolution();
}

inline const BinTable &PixelSelectorAll::bins() const noexcept { return _selectors.front().bins(); }

template <typename N>
inline bool PixelSelectorAll::iterator<N>::Pair::operator<(const Pair &other) const noexcept {
  return first < other.first;
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::Pair::operator>(const Pair &other) const noexcept {
  return first > other.first;
}

template <typename N>
inline PixelSelectorAll::iterator<N>::iterator(const PixelSelectorAll &selector)
    : _selectors(std::make_shared<SelectorQueue>()),
      _its(std::make_shared<ItPQueue>()),
      _buff(std::make_shared<std::vector<Pixel<N>>>()) {
  std::for_each(selector._selectors.begin(), selector._selectors.end(),
                [&](const PixelSelector &sel) { _selectors->push(&sel); });

  _chrom1_id = _selectors->front()->chrom1().id();
  init_iterators();
  read_next_chunk();
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::operator==(const iterator<N> &other) const noexcept {
  if (!_buff || !other._buff) {
    return _buff == other._buff;
  }

  assert(_i < _buff->size());
  assert(other._i < other._buff->size());
  return (*_buff)[_i] == (*other._buff)[_i];
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::operator!=(const iterator<N> &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator*() const -> const_reference {
  assert(_buff);
  assert(_i < _buff->size());
  return (*_buff)[_i];
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator->() const -> const_pointer {
  return &*(*this);
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator++() -> iterator & {
  assert(_buff);
  if (++_i == _buff->size()) {
    read_next_chunk();
  }
  return *this;
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N>
inline void PixelSelectorAll::iterator<N>::init_iterators() {
  assert(_its->empty());
  if (_selectors->empty()) {
    _buff = nullptr;
    return;
  }

  if (_its.use_count() != 1) {
    _its = std::make_shared<ItPQueue>();
  }

  while (!_selectors->empty() && _selectors->front()->chrom1().id() == _chrom1_id) {
    auto *sel = _selectors->front();
    _selectors->pop();
    _its->emplace(Pair{sel->begin<N>(), sel->end<N>()});
  }
}

template <typename N>
inline void PixelSelectorAll::iterator<N>::read_next_chunk() {
  if (_selectors->empty()) {
    _buff = nullptr;  // signal end
    return;
  }

  if (_its->empty()) {
    ++_chrom1_id;
    init_iterators();
    return read_next_chunk();
  }

  auto [first, last] = _its->top();
  _its->pop();

  if (first == last) {
    return read_next_chunk();
  }

  if (_buff.use_count() != 1) {
    _buff = std::make_shared<std::vector<Pixel<N>>>();
  }
  _buff->clear();
  _i = 0;

  const auto bin1 = first->coords.bin1;
  while (first != last && first->coords.bin1 == bin1) {
    _buff->push_back(*first++);
  }
  _its->emplace(Pair{first, last});
}

}  // namespace hictk::hic
