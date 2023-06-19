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

inline std::size_t PixelSelector::estimate_optimal_cache_size() const {
  if (_reader.index().empty()) {
    return 0;  // should we throw instead?
  }

  std::seed_seq sseq({_reader.index().size()});
  std::mt19937_64 rand_eng(sseq);

  // Find block with the largest compressed size and use it to guess the compression ratio
  const auto idx =
      std::max_element(_reader.index().begin(), _reader.index().end(),
                       [&](const internal::BlockIndex &idx1, const internal::BlockIndex &idx2) {
                         return idx1.compressed_size_bytes() < idx2.compressed_size_bytes();
                       });

  const auto num_pixels_in_largest_block = _reader.read(chrom1(), chrom2(), *idx)->size();
  const auto compression_ratio = static_cast<std::size_t>(
      std::ceil(double(num_pixels_in_largest_block) / double(idx->compressed_size_bytes())));

  // Try to guess the average block size
  std::size_t avg_block_size = num_pixels_in_largest_block;

  auto first_idx = _reader.index().begin();
  auto last_idx = _reader.index().end();
  auto samples = (std::min(100UL, _reader.index().size() - 1));
  // index is backed by a hashmap, so iteration should be somewhat random
  for (std::size_t i = 0; i < samples && first_idx++ != last_idx; ++i) {
    avg_block_size += first_idx->compressed_size_bytes() * compression_ratio;
  }
  avg_block_size /= samples + 1;

  // Try to guess how many blocks overlap a single row of pixels
  std::size_t avg_blocks_per_row = 0;
  const auto bin_size = bins().bin_size();

  const std::size_t first_bin_id = 0;
  const std::size_t last_bin_id =
      bins().at(coord1().bin1.chrom(), coord1().bin1.chrom().size()).rel_id() - 1;
  samples = 10;
  for (std::size_t i = 0; i < samples; ++i) {
    const auto bin_id =
        std::uniform_int_distribution<std::size_t>{first_bin_id, last_bin_id}(rand_eng);

    const auto pos1 = static_cast<std::uint32_t>(bin_id * bin_size);
    const auto pos2 = (std::min)(pos1 + bin_size, coord1().bin1.chrom().size());

    const auto coord1_ = PixelCoordinates(bins().at(coord1().bin1.chrom(), pos1),
                                          bins().at(coord1().bin1.chrom(), pos2));

    std::vector<internal::BlockIndex> buffer{};
    _reader.index().find_overlaps(coord1_, coord2(), buffer, true);
    avg_blocks_per_row += buffer.size();
  }
  avg_blocks_per_row /= samples;

  return avg_blocks_per_row * avg_block_size;
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
  // clang-format off
  return _sel == other._sel       &&
         size() == other.size();
  // clang-format on
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<(const iterator &other) const noexcept {
  assert(_sel == other._sel);
  assert(!!_sel);
  return _pixels_processed < other._pixels_processed;
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
inline const std::vector<internal::BlockIndex>
    &PixelSelector::iterator<N>::find_blocks_overlapping_next_chunk(std::size_t num_bins) {
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
inline PixelSelectorAll::iterator<N>::iterator(const PixelSelectorAll &sel) : _sel(&sel) {
  std::vector<PixelSelector::iterator<N>> heads;
  std::vector<PixelSelector::iterator<N>> tails;

  _it = _sel->_selectors.begin();
  setup_next_pixel_merger();
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::operator==(const iterator<N> &other) const noexcept {
  return _value == other._value;
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::operator!=(const iterator<N> &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator*() const -> const_reference {
  return _value;
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator->() const -> const_pointer {
  return &_value;
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator++() -> iterator & {
  _value = _merger->next();
  if (!_value && _it != _sel->_selectors.end()) {
    setup_next_pixel_merger();
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
inline void PixelSelectorAll::iterator<N>::setup_next_pixel_merger() {
  assert(_it != _sel->_selectors.end());

  if (_it != _sel->_selectors.begin()) {
    std::for_each(_sel->_selectors.begin(), _it - 1,
                  [](const auto &sel) { sel.evict_blocks_from_cache(); });
  }

  auto chrom1 = _it->chrom1();
  auto first_sel = _it;
  auto last_sel = std::find_if(first_sel, _sel->_selectors.end(),
                               [&](const PixelSelector &s) { return s.chrom1() != chrom1; });

  std::vector<PixelSelector::iterator<N>> heads;
  std::vector<PixelSelector::iterator<N>> tails;

  while (first_sel != last_sel) {
    const auto &sel = *first_sel;
    if (sel.chrom1() != chrom1) {
      if (!heads.empty()) {
        break;
      } else {
        chrom1 = sel.chrom1();
      }
    }
    auto first = sel.template begin<N>();
    auto last = sel.template end<N>();
    if (first != last) {
      heads.emplace_back(std::move(first));
      tails.emplace_back(std::move(last));
    }
    ++first_sel;
  }

  _it = last_sel;

  if (heads.empty()) {
    *this = iterator{};
    return;
  }

  _merger = std::make_shared<PixelMerger>(std::move(heads), std::move(tails));
  _value = _merger->next();

  if (!_value) {
    *this = iterator{};
    return;
  }
}

}  // namespace hictk::hic
