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
    : _reader(std::move(hfs_), footer_->index(), std::move(bins_), std::move(cache_), coord1_,
              coord2_),
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

inline SerializedPixel PixelSelector::process_interaction(SerializedPixel record) const {
  const auto &c1Norm = _footer->c1Norm();
  const auto &c2Norm = _footer->c2Norm();
  const auto &expected = _footer->expectedValues();

  assert(is_inter() || record.bin1_id <= record.bin2_id);

  const auto skipNormalization =
      normalization() == NormalizationMethod::NONE || matrix_type() == MatrixType::expected;

  if (!skipNormalization) {
    const auto bin1 = static_cast<std::size_t>(record.bin1_id);
    const auto bin2 = static_cast<std::size_t>(record.bin2_id);
    assert(bin1 < c1Norm.size());
    assert(bin2 < c2Norm.size());
    record.count /= static_cast<float>(c1Norm[bin1] * c2Norm[bin2]);
  }

  record.bin1_id *= resolution();
  record.bin2_id *= resolution();

  if (matrix_type() == MatrixType::observed) {
    return record;
  }

  const auto expectedCount = [&]() {
    if (is_inter()) {
      return float(_reader.avg());
    }

    const auto i = static_cast<std::size_t>((record.bin2_id - record.bin1_id) / resolution());
    assert(i < expected.size());
    return float(expected[i]);
  }();

  if (matrix_type() == MatrixType::expected) {
    record.count = expectedCount;
    return record;
  }

  assert(matrix_type() == MatrixType::oe);
  record.count /= expectedCount;

  return record;
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all() const {
  return {begin<N>(), end<N>()};
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all_dbg() const {
  std::vector<Pixel<N>> buffer{};

  auto bin1 = coord1().bin1.rel_id();
  auto bin2 = coord1().bin2.rel_id() + 1;
  auto bin3 = coord2().bin1.rel_id();
  auto bin4 = coord2().bin2.rel_id() + 1;

  std::size_t i = 0;

  for (const auto &block_idx : _reader.index()) {
    auto blk = _reader.read(block_idx);
    if (!blk) {
      continue;
    }
    std::ofstream ofs(fmt::format(FMT_STRING("/tmp/test{}.bed"), i));
    // Obs we use open-closed interval instead of open-open like is done in straw
    for (const auto &[b1, row] : blk->find_overlap(bin1, bin2)) {
      for (const auto &tp : row) {
        const auto &b2 = tp.bin2_id;
        auto tmp_record = process_interaction(SerializedPixel{
            static_cast<std::int64_t>(b1), static_cast<std::int64_t>(b2), tp.count});
        fmt::print(ofs, FMT_STRING("{}\t{}\t{}\t{}\n"), i, tmp_record.bin1_id, tmp_record.bin2_id,
                   blk->id());
      }

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
              PixelCoordinates{bins().at(chrom1(), static_cast<std::uint32_t>(record.bin1_id)),
                               bins().at(chrom2(), static_cast<std::uint32_t>(record.bin2_id))},
              record.count);
        }
      }
    }
    ++i;
  }
  // Only interactions from the same block are guaranteed to already be sorted
  std::sort(buffer.begin(), buffer.end());
  return buffer;
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
inline bool PixelSelector::iterator<N>::discard() const noexcept {
  if (is_at_end()) {
    return true;
  }

  assert(!!_buffer);
  if (_buffer->empty()) {
    return true;
  }

  const auto &pixel = _buffer->front();
  // clang-format off
  return pixel.coords.bin1 < coord1().bin1 ||
         pixel.coords.bin1 > coord1().bin2 ||
         pixel.coords.bin2 < coord2().bin1 ||
         pixel.coords.bin2 > coord2().bin2;
  // clang-format on
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
inline internal::Index PixelSelector::iterator<N>::find_blocks_overlapping_current_row() {
  const auto end_pos = coord2().bin2.start();
  const auto pos1 = (std::min)(end_pos, static_cast<std::uint32_t>(_bin1_id) * bins().bin_size());
  const auto pos2 = (std::min)(end_pos, pos1 + bins().bin_size());

  const auto coord1_ = PixelCoordinates(bins().at(coord1().bin1.chrom(), pos1),
                                        bins().at(coord1().bin1.chrom(), pos2));

  return _sel->_reader.index().subset(coord1_, coord2());
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
  const auto bin1 =
      bins().at(coord1().bin1.chrom(), static_cast<std::uint32_t>(_bin1_id) * bin_size);
  bool pixels_are_sorted = true;
  for (const auto block_idx : blocks) {
    const auto blk = _sel->_reader.read(block_idx);
    const auto match = blk->find(_bin1_id);
    if (match == blk->end()) {
      continue;
    }

    const auto &pixels = match->second;
    auto first = std::lower_bound(pixels.begin(), pixels.end(), coord2().bin1.rel_id(),
                                  [](const internal::InteractionBlock::ThinPixel &pixel,
                                     std::size_t bin_id) { return pixel.bin2_id < bin_id; });

    while (first != pixels.end()) {
      if (first->bin2_id > coord2().bin2.rel_id()) {
        break;
      }
      const auto pos2 = static_cast<std::uint32_t>(first->bin2_id) * bin_size;
      if constexpr (std::is_integral_v<N>) {
        _buffer->emplace_back(
            Pixel<N>{PixelCoordinates{bin1, bins().at(coord2().bin1.chrom(), pos2)},
                     static_cast<N>(std::round(first->count))});
      } else {
        _buffer->emplace_back(
            Pixel<N>{PixelCoordinates{bin1, bins().at(coord2().bin1.chrom(), pos2)},
                     conditional_static_cast<N>(first->count)});
      }

      pixels_are_sorted &= _buffer->size() < 2 || *(_buffer->end() - 2) < _buffer->back();
      ++first;
    }
  }
  if (!pixels_are_sorted) {
    std::sort(_buffer->begin(), _buffer->end());
  }
  _bin1_id++;
}

}  // namespace hictk::hic
