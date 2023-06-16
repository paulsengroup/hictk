// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
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
                                    std::shared_ptr<internal::BlockLRUCache> cache_,
                                    std::shared_ptr<const BinTable> bins_,
                                    PixelCoordinates coords) noexcept
    : PixelSelector(std::move(hfs_), std::move(footer_), std::move(cache_), std::move(bins_),
                    coords, std::move(coords)) {}

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockLRUCache> cache_,
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

inline internal::InteractionBlock::ThinPixel PixelSelector::transform_pixel(
    std::size_t bin1, internal::InteractionBlock::ThinPixel pixel) const {
  const auto &c1Norm = _footer->c1Norm();
  const auto &c2Norm = _footer->c2Norm();
  const auto &expected = _footer->expectedValues();

  assert(is_inter() || bin1 <= pixel.bin2_id);

  const auto skipNormalization =
      normalization() == NormalizationMethod::NONE || matrix_type() == MatrixType::expected;

  if (!skipNormalization) {
    const auto bin2 = pixel.bin2_id;
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

    const auto i = (pixel.bin2_id - bin1);
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

template <typename N>
inline PixelSelector::iterator<N>::iterator(const PixelSelector &sel)
    : _sel(&sel), _bin1_id(coord1().bin1.rel_id()), _buffer(std::make_shared<BufferT>()) {
  if (_sel->_reader.index().empty()) {
    *this = at_end(sel);
    return;
  }

  while (_buffer->empty()) {
    read_next_row();
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
    read_next_row();
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
inline std::vector<internal::BlockIndex>
PixelSelector::iterator<N>::find_blocks_overlapping_current_row() {
  const auto end_pos = coord1().bin2.start();
  const auto pos1 = (std::min)(end_pos, static_cast<std::uint32_t>(_bin1_id) * bins().bin_size());
  const auto pos2 = (std::min)(end_pos, pos1 + bins().bin_size());

  const auto coord1_ = PixelCoordinates(bins().at(coord1().bin1.chrom(), pos1),
                                        bins().at(coord1().bin1.chrom(), pos2));

  return _sel->_reader.index().find_overlaps(coord1_, coord2());
}

template <typename N>
inline void PixelSelector::iterator<N>::read_next_row() {
  assert(!!_sel);
  const auto blocks = find_blocks_overlapping_current_row();
  if (blocks.empty() || _bin1_id > coord1().bin2.rel_id()) {
    *this = at_end(*_sel);
    return;
  }

  if (_buffer.use_count() != 1) {
    _buffer = std::make_shared<BufferT>(_buffer->capacity());
  }

  _buffer->clear();
  _buffer_i = 0;
  const auto bin_size = bins().bin_size();
  const auto &chrom1 = coord1().bin1.chrom();
  const auto bin1 = bins().at(chrom1, static_cast<std::uint32_t>(_bin1_id) * bin_size);
  for (const auto block_idx : blocks) {
    const auto blk = _sel->_reader.read(chrom1, coord2().bin1.chrom(), block_idx);
    const auto match = blk->find(_bin1_id);
    if (match == blk->end()) {
      continue;
    }

    const auto &pixels = match->second;
    auto first = std::lower_bound(pixels.begin(), pixels.end(), coord2().bin1.rel_id(),
                                  [](const internal::InteractionBlock::ThinPixel &pixel,
                                     std::size_t bin_id) { return pixel.bin2_id < bin_id; });
    while (first != pixels.end()) {
      const auto p = _sel->transform_pixel(_bin1_id, *first);
      if (p.bin2_id > coord2().bin2.rel_id()) {
        break;
      }

      const auto pos2 = static_cast<std::uint32_t>(p.bin2_id) * bin_size;
      if constexpr (std::is_integral_v<N>) {
        _buffer->emplace_back(
            Pixel<N>{PixelCoordinates{bin1, bins().at(coord2().bin1.chrom(), pos2)},
                     conditional_static_cast<N>(std::round(p.count))});
      } else {
        _buffer->emplace_back(
            Pixel<N>{PixelCoordinates{bin1, bins().at(coord2().bin1.chrom(), pos2)},
                     conditional_static_cast<N>(p.count)});
      }
      ++first;
    }
  }
  assert(std::is_sorted(_buffer->begin(), _buffer->end()));
  _bin1_id++;
}

inline PixelSelectorAll::PixelSelectorAll(std::vector<PixelSelector> selectors_) noexcept
    : _selectors(std::move(selectors_)) {}

template <typename N>
inline auto PixelSelectorAll::begin() const -> iterator<N> {
  return cbegin<N>();
}
template <typename N>
inline auto PixelSelectorAll::cbegin() const -> iterator<N> {
  return iterator<N>(_selectors);
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
inline PixelSelectorAll::iterator<N>::iterator(const std::vector<PixelSelector> &selectors_) {
  std::vector<PixelSelector::iterator<N>> heads;
  std::vector<PixelSelector::iterator<N>> tails;

  for (const auto &sel : selectors_) {
    auto first = sel.begin<N>();
    auto last = sel.end<N>();
    if (first != last) {
      heads.emplace_back(std::move(first));
      tails.emplace_back(std::move(last));
    }
  }

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

template <typename N>
inline bool PixelSelectorAll::iterator<N>::operator==(const iterator<N> &other) const noexcept {
  return _i == other._i && _value == other._value;
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
  if (!_value) {
    *this = iterator{};
  }
  return *this;
}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

}  // namespace hictk::hic
