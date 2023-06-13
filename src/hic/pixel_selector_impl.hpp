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

  for (const auto &block_idx : _reader.grid()) {
    auto blk = _reader.read(*block_idx.block_idx);
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
              PixelCoordinates{
                  bins().at(_footer->chrom1(), static_cast<std::uint32_t>(record.bin1_id)),
                  bins().at(_footer->chrom2(), static_cast<std::uint32_t>(record.bin2_id))},
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
    : _sel(&sel),
      _grid(std::make_shared<internal::BlockGrid>(_sel->_reader.grid())),
      _idx(_grid->begin()),
      _bin1_id(coord1().bin1.rel_id()),
      _buffer(std::make_shared<BufferT>(1000)) {
  if (_grid->size() == 0) {
    *this = at_end(sel);
    return;
  }
  auto blk = read_block();
  read_chunk_of_pixels(*blk);
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
         _idx == other._idx       &&
         size() == other.size();
  // clang-format on
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
  if (_idx != other._idx) {
    return _idx->block_idx->id < other._idx._node->block_idx->id;
  }
  if (_bin1_id != other._bin1_id) {
    return _bin1_id < other._bin1_id;
  }
  return size() < other.size();
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

  if (!!_buffer) {
    ++_buffer_i;
  }

  while (_buffer_i == size()) {
    seek_to_next_block();
    if (is_at_end()) {
      break;
    }
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
inline void PixelSelector::iterator<N>::seek_to_next_block() {
  auto blk = read_block();
  if (!blk) {
    _buffer = nullptr;
    return;
  }
  assert(blk);
  const auto block_was_fully_read = _bin1_id == (--blk->end())->first;
  if (block_was_fully_read) {
    mark_block_as_fully_read();
  }

  assert(_idx < _grid->end());
  assert(_idx->current_row < _grid->end());
  const auto end_of_row = (_idx + 1 == _grid->end()) || _idx->row != (_idx + 1)->row;
  if (end_of_row) {
    if (!_grid) {
      *this = at_end(*_sel);
    }
    fmt::print(FMT_STRING("next row...\n"));
    _idx = std::find_if(_idx->current_row, _idx + 1,
                        [](const auto &node) { return node.block_idx != nullptr; });
    ++_bin1_id;
  } else {
    fmt::print(FMT_STRING("next col...\n"));
    ++_idx;
  }

  if (_idx == _grid->end()) {
    *this = at_end(*_sel);
    return;
  }

  assert(_idx->block_idx);
  blk = _sel->_reader.read(*_idx->block_idx);
  if (!blk) {
    *this = at_end(*_sel);
    return;
  }

  if (_bin1_id != blk->at(_bin1_id)->first) {
    seek_to_next_block();
  }
  read_chunk_of_pixels(*blk);
}

template <typename N>
inline void PixelSelector::iterator<N>::mark_block_as_fully_read() {
  // deal with calls to operator++(int)
  const auto i = std::distance(_grid->begin(), _idx);
  _grid = std::make_shared<internal::BlockGrid>(*_grid);
  _idx = _grid->begin() + i;
  _idx->block_idx = nullptr;
}

template <typename N>
inline std::shared_ptr<const internal::InteractionBlock>
PixelSelector::iterator<N>::read_block() noexcept {
  if (!_grid) {
    return nullptr;
  }
  assert(_sel);
  assert(_idx != _grid->end());
  assert(_idx->block_idx);
  return _sel->_reader.read(*_idx->block_idx);
}

template <typename N>
inline void PixelSelector::iterator<N>::read_chunk_of_pixels(
    const internal::InteractionBlock &blk) {
  auto row_it = blk.at(_bin1_id);
  if (_bin1_id != row_it->first) {
    return;
  }

  auto first = row_it->second.begin();
  auto last = row_it->second.end();

  first =
      std::lower_bound(first, last, coord2().bin1.start(),
                       [&](const internal::InteractionBlock::ThinPixel &pixel, const auto &pos) {
                         return pixel.bin2_id * bins().bin_size() < pos;
                       });
  // clang-format off
  const auto buff_capacity =
      std::min({!!_buffer ? _buffer->capacity() : std::size_t(1000),
                static_cast<std::size_t>(std::distance(first, last))});
  // clang-format on
  // This is fine, as iterators are not thread safe anyway
  if (_buffer.use_count() == 1) {
    _buffer->resize(buff_capacity);
  } else {
    _buffer = std::make_shared<BufferT>(buff_capacity);
  }
  _buffer->clear();
  _buffer_i = 0;

  const auto pos1 = static_cast<std::uint32_t>(_bin1_id) * bins().bin_size();
  const auto bin1 = bins().at(coord1().bin1.chrom(), pos1);
  do {
    const auto pos2 = static_cast<std::uint32_t>(first->bin2_id) * bins().bin_size();
    auto bin2 = bins().at(coord2().bin1.chrom(), pos2);

    if (bin2 > coord2().bin2) {
      break;
    }

    if constexpr (std::is_integral_v<N>) {
      const auto count = static_cast<N>(std::round(first->count));
      _buffer->emplace_back(bin1, std::move(bin2), count);
    } else {
      const auto count = conditional_static_cast<N>(first->count);
      _buffer->emplace_back(bin1, std::move(bin2), count);
    }
  } while (++first != last);

  if (first != last && _bin1_id + 1 > coord1().bin2.rel_id()) {
    _grid = nullptr;
    _idx = {};
  }
}
}  // namespace hictk::hic
