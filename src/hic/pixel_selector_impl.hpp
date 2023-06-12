// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "hictk/common.hpp"
#include "hictk/fmt.hpp"  // TODO: remove me
#include "hictk/hic/hic_file_stream.hpp"

namespace hictk::hic {

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileStream> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockLRUCache> cache_,
                                    std::shared_ptr<const BinTable> bins_,
                                    PixelCoordinates coords) noexcept
    : PixelSelector(std::move(hfs_), std::move(footer_), std::move(cache_), std::move(bins_),
                    coords, std::move(coords)) {}

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileStream> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockLRUCache> cache_,
                                    std::shared_ptr<const BinTable> bins_, PixelCoordinates coord1_,
                                    PixelCoordinates coord2_) noexcept
    : _reader(std::move(hfs_), *footer_, std::move(bins_), std::move(cache_), coord1_, coord2_),
      _footer(std::move(footer_)),
      _coord1(std::move(coord1_)),
      _coord2(std::move(coord2_)) {}

inline bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  return _footer == other._footer && _coord1 == other._coord1 && _coord2 == other._coord2;
}

inline bool PixelSelector::operator!=(const PixelSelector &other) const noexcept {
  return !(*this == other);
}
/*
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
*/

inline SerializedPixel PixelSelector::process_interaction(SerializedPixel record) const {
  const auto &c1Norm = _footer->c1Norm();
  const auto &c2Norm = _footer->c2Norm();
  const auto &expected = _footer->expectedValues();

  assert(is_inter() || record.bin1_id <= record.bin2_id);

  const auto skipNormalization = _footer->normalization() == NormalizationMethod::NONE ||
                                 _footer->matrix_type() == MatrixType::expected;

  if (!skipNormalization) {
    const auto bin1 = static_cast<std::size_t>(record.bin1_id);
    const auto bin2 = static_cast<std::size_t>(record.bin2_id);
    assert(bin1 < c1Norm.size());
    assert(bin2 < c2Norm.size());
    record.count /= static_cast<float>(c1Norm[bin1] * c2Norm[bin2]);
  }

  record.bin1_id *= _footer->resolution();
  record.bin2_id *= _footer->resolution();

  if (_footer->matrix_type() == MatrixType::observed) {
    return record;
  }

  const auto expectedCount = [&]() {
    if (is_inter()) {
      return float(_reader.avg());
    }

    const auto i =
        static_cast<std::size_t>((record.bin2_id - record.bin1_id) / _footer->resolution());
    assert(i < expected.size());
    return float(expected[i]);
  }();

  if (_footer->matrix_type() == MatrixType::expected) {
    record.count = expectedCount;
    return record;
  }

  assert(_footer->matrix_type() == MatrixType::oe);
  record.count /= expectedCount;

  return record;
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all() const {
  // return {begin<N>(), end<N>()};
  std::vector<Pixel<N>> buffer{};
  std::size_t empty_blocks = 0;

  auto bin1 = coord1().bin1.rel_id();
  auto bin2 = coord1().bin2.rel_id() + 1;
  auto bin3 = coord2().bin1.rel_id();
  auto bin4 = coord2().bin2.rel_id() + 1;

  for (const auto &block_idx : _reader.grid()) {
    auto blk = _reader.read(*block_idx.block);
    if (!blk) {
      empty_blocks++;
      continue;
    }

    // Obs we use open-closed interval instead of open-open like is done in straw
    for (const auto &[b1, row] : blk->find_overlap(bin1, bin2)) {
      if (b1 >= bin2) {
        // We're past the last row overlapping the query
        break;
      }
      for (const auto &tp : row) {
        const auto &b2 = tp.bin2_id;
        if (b1 < bin1 || b2 < bin3) {
          // We're upstream of the first column overlapping the query (if any)
          continue;
        }

        if (b2 >= bin4) {
          // We're past the last column overlapping the query for the current row
          break;
        }

        auto record = process_interaction(SerializedPixel{static_cast<std::int64_t>(b1),
                                                          static_cast<std::int64_t>(b2), tp.count});
        if (std::isfinite(record.count)) {
          buffer.emplace_back(
              PixelCoordinates{
                  bins().at(_footer->chrom1(), static_cast<std::uint32_t>(record.bin1_id)),
                  bins().at(_footer->chrom2(), static_cast<std::uint32_t>(record.bin2_id))},
              record.count);
        }
      }
    }
  }
  // Only interactions from the same block are guaranteed to already be sorted
  std::sort(buffer.begin(), buffer.end());
  return buffer;
}

inline const PixelCoordinates &PixelSelector::coord1() const noexcept { return _coord1; }
inline const PixelCoordinates &PixelSelector::coord2() const noexcept { return _coord2; }

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
/*
template <typename N>
inline PixelSelector::iterator<N>::iterator(const PixelSelector &sel)
    : _sel(&sel), _blk(_sel->_reader.read(this->coord1(), coord2())) {
  if (!_blk) {
    *this = at_end(sel);
    return;
  }

  _row = _blk->begin();
  seek_to_next_overlap();
  _value = *(*this);
}

template <typename N>
inline auto PixelSelector::iterator<N>::at_end(const PixelSelector &sel) -> iterator<N> {
  iterator it{};

  it._sel = &sel;
  it._blk = nullptr;  // This signals we're at end

  return it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator==(const iterator &other) const noexcept {
  return _sel == other._sel && _blk == other._blk && _row == other._row && _col_i == other._col_i;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<(const iterator &other) const noexcept {
  assert(!!_sel);
  assert(_sel->coord1() == other._sel->coord1());
  assert(_sel->coord2() == other._sel->coord2());
  return _value < other._value;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator<=(const iterator &other) const noexcept {
  assert(!!_sel);
  assert(_sel->coord1() == other._sel->coord1());
  assert(_sel->coord2() == other._sel->coord2());
  return _value <= other._value;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator>(const iterator &other) const noexcept {
  assert(!!_sel);
  assert(_sel->coord1() == other._sel->coord1());
  assert(_sel->coord2() == other._sel->coord2());
  return _value > other._value;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator>=(const iterator &other) const noexcept {
  assert(!!_sel);
  assert(_sel->coord1() == other._sel->coord1());
  assert(_sel->coord2() == other._sel->coord2());
  return _value >= other._value;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator*() const -> const_reference {
  assert(!!_sel);

  _value.coords = {bins().at(coord1().bin1.chrom(), pos1()),
                   bins().at(coord2().bin1.chrom(), pos2())};
  _value.count = count();
  return _value;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++() -> iterator & {
  assert(!!_sel);
  assert(!is_at_end());

  if (discard()) {
    this->_col_i--;
    fmt::print(FMT_STRING("operator++={}\n"), *(*this));
    this->_col_i++;
    seek_to_next_block();
  } else {
    fmt::print(FMT_STRING("operator++{}\n"), *(*this));
    assert(_col_i < row().size());
    _col_i++;
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
inline bool PixelSelector::iterator<N>::discard() const noexcept {
  if (is_at_end()) {
    return false;
  }
  // clang-format off
  return _col_i >= row().size() ||
         pos1() < coord1().bin1.start() ||
         pos1() > coord1().bin2.start() ||
         pos2() < coord2().bin1.start() ||
         pos2() > coord2().bin2.start();
  // clang-format on
}

template <typename N>
inline bool PixelSelector::iterator<N>::is_at_end() const noexcept {
  return !_blk;
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
inline const internal::InteractionBlock::Row &PixelSelector::iterator<N>::row() const noexcept {
  assert(!is_at_end());
  return _row->second;
}

template <typename N>
inline std::uint32_t PixelSelector::iterator<N>::pos1() const noexcept {
  assert(!is_at_end());
  return static_cast<std::uint32_t>(_row->first) * bins().bin_size();
}

template <typename N>
inline std::uint32_t PixelSelector::iterator<N>::pos2() const noexcept {
  assert(!is_at_end());
  assert(_col_i < row().size());
  const auto &thin_pixel = row()[_col_i];
  return static_cast<std::uint32_t>(thin_pixel.bin2_id) * bins().bin_size();
}

template <typename N>
inline N PixelSelector::iterator<N>::count() const noexcept {
  assert(!is_at_end());
  assert(_col_i < row().size());
  const auto &thin_pixel = row()[_col_i];
  if constexpr (std::is_integral_v<N>) {
    return static_cast<N>(std::round(thin_pixel.count));
  } else {
    return conditional_static_cast<N>(thin_pixel.count);
  }
}

template <typename N>
inline void PixelSelector::iterator<N>::seek_to_next_block() {
  assert(!is_at_end());

  // Figure out whether we should move right or down

  if (_value.coords.bin1.rel_id() == (--_blk->end())->first) {
    // Advance row and reset col
    // TODO
    assert(false);
    // const auto pos1 = (static_cast<std::uint32_t>(_blk->last_row() + 2) * bins().bin_size()) -
    // 1; const PixelCoordinates coords_{bins().at(coord1().bin1.chrom(), pos1), coord1().bin2};
    // auto blk = _sel->_reader.read(coords_, coord2());
    // assert(blk != _blk);
    // _blk = blk;
  } else {
    // Advance col and leave row alone
    _blk = _sel->_reader.read_after(_blk->id());
  }

  // fmt::print(FMT_STRING("jumping to {}:{}\n"), coords1_, coord2());
  if (!_blk) {
    *this = at_end(*_sel);
    return;
  }
  _col_i = 0;
  _row = _blk->begin();

  seek_to_next_overlap();
  _value = *(*this);
}

template <typename N>
inline void PixelSelector::iterator<N>::seek_to_next_overlap() noexcept {
  assert(!is_at_end());

  for (; _row != _blk->end(); ++_row) {
    if (pos1() >= coord1().bin1.start()) {
      break;
    }
  }

  for (_col_i = 0; _col_i < row().size(); ++_col_i) {
    if (pos2() >= coord2().bin1.start()) {
      return;
    }
  }
}
 */

}  // namespace hictk::hic
