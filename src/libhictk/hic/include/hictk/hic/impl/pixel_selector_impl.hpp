// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <memory>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic {

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockCache> cache_,
                                    std::shared_ptr<const BinTable> bins_,
                                    const PixelCoordinates &coords,
                                    std::optional<std::uint64_t> diagonal_band_width)
    : PixelSelector(std::move(hfs_), std::move(footer_), std::move(cache_), std::move(bins_),
                    coords, coords, diagonal_band_width) {}

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<internal::BlockCache> cache_,
                                    std::shared_ptr<const BinTable> bins_, PixelCoordinates coord1_,
                                    PixelCoordinates coord2_,
                                    std::optional<std::uint64_t> diagonal_band_width)
    : _reader(std::make_shared<internal::HiCBlockReader>(std::move(hfs_), footer_->index(),
                                                         std::move(bins_), std::move(cache_))),
      _footer(std::move(footer_)),
      _coord1(std::make_shared<const PixelCoordinates>(std::move(coord1_))),
      _coord2(std::make_shared<const PixelCoordinates>(std::move(coord2_))),
      _diagonal_band_width(diagonal_band_width) {
  const auto query_is_cis = _coord1->bin1.chrom() == _coord2->bin1.chrom();
  if ((!query_is_cis && _coord1->bin1 > _coord2->bin1) ||
      (query_is_cis && _coord1->bin1.start() > _coord2->bin1.start())) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("query {}:{}-{}; {}:{}-{}; overlaps with the lower-triangle of the matrix"),
        _coord1->bin1.chrom().name(), _coord1->bin1.start(), _coord1->bin2.end(),
        _coord2->bin1.chrom().name(), _coord2->bin1.start(), _coord2->bin2.end()));
  }
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline PixelSelector::~PixelSelector() noexcept {
  try {
    if (_reader) {
      clear_cache();
    }
  } catch (const std::exception &e) {
    SPDLOG_WARN(
        FMT_STRING("failed to clear the PixelSelector interaction cache of file \"{}\". "
                   "Applications are not affected. However, if this happens repeatedly the system "
                   "may run out of available memory. Reason: {}"),
        _footer->path(), e.what());
  } catch (...) {
    SPDLOG_WARN(
        FMT_STRING("failed to clear the PixelSelector interaction cache of file \"{}\". "
                   "Applications are not affected. However, if this happens repeatedly the system "
                   "may run out of available memory."),
        _footer->path());
  }
}

inline bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  return _reader->index().chrom1() == _reader->index().chrom2() && *_coord1 == *other._coord1 &&
         *_coord2 == *other._coord2;
}

inline bool PixelSelector::operator!=(const PixelSelector &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelector::cbegin(bool sorted) const -> iterator<N> {
  return iterator<N>(*this, sorted, _diagonal_band_width);
}

template <typename N>
inline auto PixelSelector::cend() const -> iterator<N> {
  return iterator<N>::at_end(_reader, _coord1, _coord2);
}

template <typename N>
inline auto PixelSelector::begin(bool sorted) const -> iterator<N> {
  return cbegin<N>(sorted);
}

template <typename N>
inline auto PixelSelector::end() const -> iterator<N> {
  return cend<N>();
}

template <typename N>
inline ThinPixel<N> PixelSelector::transform_pixel(ThinPixel<float> pixel) const {
  auto return_pixel = [&]() -> ThinPixel<N> {
    if constexpr (std::is_floating_point_v<N>) {
      return {pixel.bin1_id, pixel.bin2_id, conditional_static_cast<N>(pixel.count)};
    } else {
      return {pixel.bin1_id, pixel.bin2_id, static_cast<N>(std::round(pixel.count))};
    }
  };

  const auto &weights1 = _footer->weights1()(balancing::Weights::Type::DIVISIVE);
  const auto &weights2 = _footer->weights2()(balancing::Weights::Type::DIVISIVE);
  const auto &expected = _footer->expectedValues();

  const auto bin1 = pixel.bin1_id;
  const auto bin2 = pixel.bin2_id;

  assert(is_inter() || bin1 <= bin2);

  const auto skipNormalization =
      normalization() == balancing::Method::NONE() || matrix_type() == MatrixType::expected;

  if (!skipNormalization) {
    assert(bin1 < weights1.size());
    assert(bin2 < weights2.size());
    pixel.count /= static_cast<float>(weights1[bin1] * weights2[bin2]);
  }

  if (matrix_type() == MatrixType::observed) {
    return return_pixel();
  }

  const auto expectedCount = [&]() {
    if (is_inter()) {
      return float(_reader->avg());
    }

    const auto i =
        std::min(bin2 - bin1, conditional_static_cast<std::uint64_t>(expected.size() - 1));
    return float(expected[i]);
  }();

  if (matrix_type() == MatrixType::expected) {
    pixel.count = expectedCount;
    return return_pixel();
  }

  assert(matrix_type() == MatrixType::oe);
  pixel.count /= expectedCount;

  return return_pixel();
}

inline PixelSelector::PixelSelector(std::shared_ptr<internal::HiCBlockReader> reader_,
                                    std::shared_ptr<const internal::HiCFooter> footer_,
                                    std::shared_ptr<const PixelCoordinates> coord1_,
                                    std::shared_ptr<const PixelCoordinates> coord2_,
                                    std::optional<std::uint64_t> diagonal_band_width)
    : _reader(std::move(reader_)),
      _footer(std::move(footer_)),
      _coord1(std::move(coord1_)),
      _coord2(std::move(coord2_)),
      _diagonal_band_width(diagonal_band_width) {}

inline PixelSelector PixelSelector::fetch(PixelCoordinates coord1_,
                                          PixelCoordinates coord2_) const {
  return {_reader, _footer, std::make_shared<const PixelCoordinates>(std::move(coord1_)),
          std::make_shared<const PixelCoordinates>(std::move(coord2_)), _diagonal_band_width};
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all() const {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel<N>> buff{};
  std::transform(begin<N>(), end<N>(), std::back_inserter(buff), [&](const ThinPixel<N> &p) {
    return Pixel<N>{{bins().at_hint(p.bin1_id, coord1().bin1.chrom()),
                     bins().at_hint(p.bin2_id, coord2().bin1.chrom())},
                    p.count};
  });
  return buff;
}

inline const PixelCoordinates &PixelSelector::coord1() const noexcept { return *_coord1; }
inline const PixelCoordinates &PixelSelector::coord2() const noexcept { return *_coord2; }
inline std::uint64_t PixelSelector::size(bool upper_triangle) const {
  if (!_coord1) {
    assert(!_coord2);
    return 0;
  }

  assert(bins().type() == BinTable::Type::fixed);

  if (_coord1->bin1.chrom() != _coord2->bin1.chrom()) {
    const auto height = _coord1->bin2.id() - _coord1->bin1.id() + 1;
    const auto width = _coord2->bin2.id() - _coord2->bin1.id() + 1;
    return height * width;
  }

  const auto start1 = _coord1->bin1.start();
  const auto start2 = _coord2->bin1.start();

  const auto end1 = (_coord1->bin2.rel_id() + 1) * resolution() + 1;
  const auto end2 = (_coord2->bin2.rel_id() + 1) * resolution() + 1;

  return area(start1, end1, start2, end2, resolution(), upper_triangle);
}

inline MatrixType PixelSelector::matrix_type() const noexcept { return metadata().matrix_type; }
inline const balancing::Method &PixelSelector::normalization() const noexcept {
  return metadata().normalization;
}
inline MatrixUnit PixelSelector::unit() const noexcept { return _reader->index().unit(); }
inline std::uint32_t PixelSelector::resolution() const noexcept {
  return _reader->index().resolution();
}

inline const Chromosome &PixelSelector::chrom1() const noexcept { return _coord1->bin1.chrom(); }
inline const Chromosome &PixelSelector::chrom2() const noexcept { return _coord2->bin1.chrom(); }

inline const balancing::Weights &PixelSelector::weights1() const noexcept {
  return _footer->weights1();
}
inline const balancing::Weights &PixelSelector::weights2() const noexcept {
  return _footer->weights2();
}

inline const BinTable &PixelSelector::bins() const noexcept { return _reader->bins(); }
inline std::shared_ptr<const BinTable> PixelSelector::bins_ptr() const noexcept {
  return _reader->bins_ptr();
}

inline const internal::HiCFooterMetadata &PixelSelector::metadata() const noexcept {
  assert(!!_footer);
  return _footer->metadata();
}

inline bool PixelSelector::is_inter() const noexcept { return !is_intra(); }

inline bool PixelSelector::is_intra() const noexcept { return chrom1() == chrom2(); }

inline bool PixelSelector::empty() const noexcept { return _reader->index().empty(); }

inline std::size_t PixelSelector::estimate_optimal_cache_size(
    [[maybe_unused]] std::size_t num_samples) const {
  if (_reader->index().empty()) {
    return 0;  // should we throw instead?
  }

  std::seed_seq sseq({_reader->index().size()});
  std::mt19937_64 rand_eng(sseq);

  // Try to guess the average block size
  std::size_t avg_block_size = 0;
  const auto &idx = _reader->index();

  std::vector<internal::BlockIndex> blocks(std::min(idx.size(), num_samples));
  std::sample(idx.begin(), idx.end(), blocks.begin(), blocks.size(), rand_eng);
  for (const auto &blki : blocks) {
    avg_block_size += _reader->read_size(chrom1(), chrom2(), blki);
  }
  avg_block_size /= blocks.size();

  // Try to guess how many blocks overlap a single row of pixels
  std::size_t max_blocks_per_row = 0;
  const auto &chrom = coord1().bin1.chrom();
  const auto bin_size = bins().resolution();

  const std::size_t first_bin_id = 0;
  const std::size_t last_bin_id = bins().at(chrom, chrom.size() - 1).rel_id();

  const auto samples = std::min(num_samples, bins().subset(chrom).size());
  for (std::size_t i = 0; i < samples; ++i) {
    const auto bin_id = std::uniform_int_distribution<std::size_t>{
        first_bin_id, std::min(last_bin_id, last_bin_id - 1)}(rand_eng);

    const auto pos1 = static_cast<std::uint32_t>(bin_id * bin_size);
    const auto bin1 = bins().at(chrom, pos1);

    auto overlap = idx.find_overlaps({bin1, bin1}, coord2(), _diagonal_band_width);
    const auto num_blocks = static_cast<std::size_t>(std::distance(overlap.begin(), overlap.end()));
    max_blocks_per_row = (std::max)(max_blocks_per_row, num_blocks);
  }

  return max_blocks_per_row * avg_block_size * sizeof(ThinPixel<float>);
}

inline void PixelSelector::clear_cache() const {
  if (_reader->cache_size() == 0) {
    return;
  }
  for (auto blki : _reader->index().find_overlaps(coord1(), coord2(), _diagonal_band_width)) {
    _reader->evict(coord1().bin1.chrom(), coord2().bin1.chrom(), blki);
  }
}

template <typename N>
inline PixelSelector::iterator<N>::iterator(const PixelSelector &sel, bool sorted,
                                            std::optional<std::uint64_t> diagonal_band_width)
    : _reader(sel._reader),
      _coord1(sel._coord1),
      _coord2(sel._coord2),
      _footer(sel._footer),
      _block_idx(preload_block_index(sel, diagonal_band_width, sorted)),
      _block_blacklist(std::make_shared<BlockBlacklist>()),
      _block_it(!!_block_idx ? _block_idx->begin() : internal::Index::Overlap::const_iterator{}),
      _buffer(std::make_shared<BufferT>()),
      _bin1_id(coord1().bin1.rel_id()),
      _diagonal_band_width(diagonal_band_width.value_or(std::numeric_limits<std::uint64_t>::max())),
      _sorted(sorted) {
  if (!!_block_idx && _block_idx->empty()) {
    *this = at_end(_reader, _coord1, _coord2);
    return;
  }

  while (!!_buffer && _buffer->empty()) {
    read_next_chunk();
  }
}

template <typename N>
inline auto PixelSelector::iterator<N>::at_end(std::shared_ptr<internal::HiCBlockReader> reader,
                                               std::shared_ptr<const PixelCoordinates> coord1,
                                               std::shared_ptr<const PixelCoordinates> coord2)
    -> iterator {
  iterator it{};

  it._reader = std::move(reader);
  it._coord1 = std::move(coord1);
  it._coord2 = std::move(coord2);
  it._buffer = nullptr;  // end of queue

  return it;
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator==(const iterator &other) const noexcept {
  const auto same_query = coord1() == other.coord1() && coord2() == other.coord2();

  if (!same_query) {
    return false;
  }

  assert(_reader->index().chrom1() == other._reader->index().chrom1());
  assert(_reader->index().chrom2() == other._reader->index().chrom2());

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
  assert(!!_reader);
  return _reader->bins();
}
template <typename N>
inline const PixelCoordinates &PixelSelector::iterator<N>::coord1() const noexcept {
  assert(!!_coord1);
  return *_coord1;
}
template <typename N>
inline const PixelCoordinates &PixelSelector::iterator<N>::coord2() const noexcept {
  assert(!!_coord2);
  return *_coord2;
}

template <typename N>
inline std::size_t PixelSelector::iterator<N>::size() const noexcept {
  return !_buffer ? 0 : _buffer->size();
}

template <typename N>
inline std::uint64_t PixelSelector::iterator<N>::bin1_id() const noexcept {
  return !is_at_end() ? (*this)->bin1_id : std::numeric_limits<std::size_t>::max();
}
template <typename N>
inline std::uint64_t PixelSelector::iterator<N>::bin2_id() const noexcept {
  return !is_at_end() ? (*this)->bin2_id : std::numeric_limits<std::size_t>::max();
}

template <typename N>
inline std::vector<internal::BlockIndex>
PixelSelector::iterator<N>::find_blocks_overlapping_next_chunk(std::size_t num_bins) const {
  const auto bin_size = bins().resolution();

  const auto end_pos = coord1().bin2.start();
  const auto pos1 = (std::min)(end_pos, _bin1_id * bins().resolution());
  const auto pos2 = (std::min)(end_pos, pos1 + static_cast<std::uint32_t>((num_bins * bin_size)));

  const auto coord1_ = PixelCoordinates(bins().at(coord1().bin1.chrom(), pos1),
                                        bins().at(coord1().bin1.chrom(), pos2));

  return _reader->index().find_overlaps(coord1_, coord2(), _diagonal_band_width);
}

template <typename N>
inline std::uint32_t PixelSelector::iterator<N>::compute_chunk_size() const noexcept {
  assert(!!_reader);
  const auto bin_size = bins().resolution();
  const auto num_bins =
      std::min(static_cast<std::uint32_t>(_reader->index().block_bin_count()),
               static_cast<std::uint32_t>(
                   0.005 * static_cast<double>(coord1().bin2.rel_id() - coord1().bin1.rel_id())));

  const auto end_pos = coord1().bin2.start();
  const auto pos1 = (std::min)(end_pos, static_cast<std::uint32_t>(_bin1_id) * bins().resolution());
  const auto pos2 = (std::min)(end_pos, pos1 + (num_bins * bin_size));

  return static_cast<std::uint32_t>((pos2 - pos1 + bin_size - 1) / bin_size);
}

template <typename N>
inline void PixelSelector::iterator<N>::read_next_chunk() {
  assert(!!_reader);
  if (!_sorted) {
    return read_next_chunk_unsorted();
  }

  const auto is_intra = coord1().bin1.chrom() == coord2().bin1.chrom();

  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  if (_reader->index().version() > 8 && is_intra) {
    return read_next_chunk_v9_intra_sorted();
  }
  read_next_chunk_sorted();
}

template <typename N>
inline void PixelSelector::iterator<N>::read_next_chunk_unsorted() {
  assert(!!_reader);

  if (_block_it == _block_idx->end()) {
    *this = at_end(_reader, _coord1, _coord2);
    return;
  }

  if (_buffer.use_count() != 1) {
    _buffer = std::make_shared<BufferT>(_buffer->capacity());
  }
  _buffer->clear();
  _buffer_i = 0;

  const auto bin1_lb = coord1().bin1.rel_id();
  const auto bin1_ub = coord1().bin2.rel_id();
  const auto bin2_lb = coord2().bin1.rel_id();
  const auto bin2_ub = coord2().bin2.rel_id();

  const auto bin1_offset = bins().at(coord1().bin1.chrom()).id();
  const auto bin2_offset = bins().at(coord2().bin1.chrom()).id();
  auto blk = *_reader->read(coord1().bin1.chrom(), coord2().bin1.chrom(), *_block_it++, false);
  for (auto p : blk) {
    if (p.bin1_id < bin1_lb || p.bin1_id > bin1_ub || p.bin2_id < bin2_lb || p.bin2_id > bin2_ub) {
      continue;
    }

    auto pt = transform_pixel(p);
    pt.bin1_id += bin1_offset;
    pt.bin2_id += bin2_offset;

    if (pt.bin2_id - pt.bin1_id >= _diagonal_band_width) {
      continue;
    }

    _buffer->emplace_back(std::move(pt));
  }
}

template <typename N>
inline void PixelSelector::iterator<N>::read_next_chunk_sorted() {
  assert(!!_reader);
  assert(_sorted);

  if (_block_it == _block_idx->end()) {
    *this = at_end(_reader, _coord1, _coord2);
    return;
  }

  if (_buffer.use_count() != 1) {
    _buffer = std::make_shared<BufferT>(_buffer->capacity());
  }
  _buffer->clear();
  _buffer_i = 0;

  const auto bin1_lb = _bin1_id;
  const auto bin1_ub = coord1().bin2.rel_id();
  const auto bin2_lb = coord2().bin1.rel_id();
  const auto bin2_ub = coord2().bin2.rel_id();

  const auto bin1_offset = bins().at(coord1().bin1.chrom()).id();
  const auto bin2_offset = bins().at(coord2().bin1.chrom()).id();
  auto first_blki = *_block_it;
  while (_block_it->coords().i1 == first_blki.coords().i1) {
    auto blk = *_reader->read(coord1().bin1.chrom(), coord2().bin1.chrom(), *_block_it, false);
    for (auto p : blk) {
      if (static_cast<std::size_t>(p.bin1_id) < bin1_lb ||
          static_cast<std::size_t>(p.bin1_id) > bin1_ub ||
          static_cast<std::size_t>(p.bin2_id) < bin2_lb ||
          static_cast<std::size_t>(p.bin2_id) > bin2_ub) {
        continue;
      }

      auto pt = transform_pixel(p);
      pt.bin1_id += bin1_offset;
      pt.bin2_id += bin2_offset;

      if (pt.bin2_id - pt.bin1_id >= _diagonal_band_width) {
        continue;
      }

      _buffer->emplace_back(std::move(pt));
    }
    if (++_block_it == _block_idx->end()) {
      break;
    }
  }

  std::sort(_buffer->begin(), _buffer->end());
}

template <typename N>
inline void PixelSelector::iterator<N>::read_next_chunk_v9_intra_sorted() {
  assert(!!_reader);
  assert(_sorted);

  if (_bin1_id > coord1().bin2.rel_id()) {
    *this = at_end(_reader, _coord1, _coord2);
    return;
  }

  if (_buffer.use_count() != 1) {
    _buffer = std::make_shared<BufferT>(_buffer->capacity());
  }

  if (_block_blacklist.use_count() != 1) {
    _block_blacklist = std::make_shared<BlockBlacklist>(*_block_blacklist);
  }
  _buffer->clear();
  _buffer_i = 0;

  const auto chunk_size = compute_chunk_size();
  const auto bin1_id_last = _bin1_id + chunk_size;

  const auto block_indexes = find_blocks_overlapping_next_chunk(chunk_size);
  for (const auto &blki : block_indexes) {
    if (_block_blacklist->contains(blki)) {
      continue;
    }

    const auto bin1_lb = _bin1_id;
    const auto bin1_ub = coord1().bin2.rel_id();
    const auto bin2_lb = coord2().bin1.rel_id();
    const auto bin2_ub = coord2().bin2.rel_id();

    const auto bin1_offset = bins().at(coord1().bin1.chrom()).id();
    const auto bin2_offset = bins().at(coord2().bin1.chrom()).id();
    bool block_overlaps_query = false;
    for (auto p : *_reader->read(coord1().bin1.chrom(), coord2().bin1.chrom(), blki)) {
      // Using bitwise operators gives a ~5% perf improvement on my machine
      // NOLINTBEGIN(*-signed-bitwise)
      using flag = std::uint_fast8_t;
      const auto pixel_overlaps_query =
          static_cast<bool>(flag(p.bin1_id >= bin1_lb) & flag(p.bin1_id <= bin1_ub) &
                            flag(p.bin2_id >= bin2_lb) & flag(p.bin2_id <= bin2_ub));
      // NOLINTEND(*-signed-bitwise)

      const auto pixel_overlaps_chunk = pixel_overlaps_query && p.bin1_id <= bin1_id_last;

      block_overlaps_query |= pixel_overlaps_query;
      if (!pixel_overlaps_chunk) {
        continue;
      }

      auto pt = transform_pixel(p);
      pt.bin1_id += bin1_offset;
      pt.bin2_id += bin2_offset;

      if (pt.bin2_id - pt.bin1_id >= _diagonal_band_width) {
        continue;
      }

      _buffer->emplace_back(std::move(pt));
    }
    if (!block_overlaps_query) {
      _reader->evict(coord1().bin1.chrom(), coord2().bin1.chrom(), blki);
      _block_blacklist->emplace(blki);
    }
  }

  std::sort(_buffer->begin(), _buffer->end());
  _bin1_id = bin1_id_last + 1;
}

template <typename N>
inline ThinPixel<N> PixelSelector::iterator<N>::transform_pixel(ThinPixel<float> pixel) const {
  assert(!!_footer);

  auto return_pixel = [&]() -> ThinPixel<N> {
    if constexpr (std::is_floating_point_v<N>) {
      return {pixel.bin1_id, pixel.bin2_id, conditional_static_cast<N>(pixel.count)};
    } else {
      return {pixel.bin1_id, pixel.bin2_id, static_cast<N>(std::round(pixel.count))};
    }
  };

  assert(_footer->weights1().type() == balancing::Weights::Type::DIVISIVE);
  assert(_footer->weights2().type() == balancing::Weights::Type::DIVISIVE);

  const auto &weights1 = _footer->weights1();
  const auto &weights2 = _footer->weights2();
  const auto &expected = _footer->expectedValues();

  const auto bin1 = pixel.bin1_id;
  const auto bin2 = pixel.bin2_id;

  const auto is_inter = coord1().bin1.chrom() != coord2().bin1.chrom();
  const auto matrix_type = _footer->matrix_type();

  assert(is_inter || bin1 <= bin2);

  const auto skip_normalization =
      _footer->normalization() == balancing::Method::NONE() || matrix_type == MatrixType::expected;

  if (!skip_normalization) {
    assert(bin1 < weights1.size());
    assert(bin2 < weights2.size());
    pixel.count /= static_cast<float>(weights1[bin1] * weights2[bin2]);
  }

  if (matrix_type == MatrixType::observed) {
    return return_pixel();
  }

  const auto expected_count = [&]() {
    if (is_inter) {
      return float(_reader->avg());
    }

    const auto i =
        std::min(bin2 - bin1, conditional_static_cast<std::uint64_t>(expected.size() - 1));
    return float(expected[i]);
  }();

  if (matrix_type == MatrixType::expected) {
    pixel.count = expected_count;
    return return_pixel();
  }

  assert(matrix_type == MatrixType::oe);
  pixel.count /= expected_count;

  return return_pixel();
}

template <typename N>
inline std::shared_ptr<const internal::Index::Overlap>
PixelSelector::iterator<N>::preload_block_index(const PixelSelector &sel,
                                                std::optional<std::uint64_t> diagonal_band_width,
                                                bool sorted) {
  const auto &idx = sel._reader->index();
  const auto is_intra = sel.coord1().bin1.chrom() == sel.coord2().bin2.chrom();
  if (sorted && is_intra && idx.version() > 8) {  // NOLINT(*-avoid-magic-numbers)
    return nullptr;
  }

  return std::make_shared<const internal::Index::Overlap>(
      idx.find_overlaps(sel.coord1(), sel.coord2(), diagonal_band_width));
}

inline PixelSelectorAll::PixelSelectorAll(
    std::vector<PixelSelector> selectors_,
    std::shared_ptr<internal::WeightCache> weight_cache) noexcept
    : _selectors(std::move(selectors_)),
      _bins(_selectors.empty() ? nullptr : _selectors.front().bins_ptr()),
      _weight_cache(std::move(weight_cache)) {}

inline bool PixelSelectorAll::empty() const noexcept { return begin<float>() == end<float>(); }

template <typename N>
inline auto PixelSelectorAll::begin(bool sorted) const -> iterator<N> {
  return cbegin<N>(sorted);
}
template <typename N>
inline auto PixelSelectorAll::cbegin(bool sorted) const -> iterator<N> {
  return iterator<N>(*this, sorted);
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
  std::transform(begin<N>(), end<N>(), std::back_inserter(buff), [&](const ThinPixel<N> &p) {
    return Pixel<N>{{bins().at(p.bin1_id), bins().at(p.bin2_id)}, p.count};
  });

  return buff;
}

inline std::uint64_t PixelSelectorAll::size(bool upper_triangle) const noexcept {
  const auto n = bins().size();
  if (upper_triangle) {
    return (n * (n + 1)) / 2;
  }
  return n * n;
}

inline MatrixType PixelSelectorAll::matrix_type() const noexcept {
  return _selectors.front().matrix_type();
}
inline const balancing::Method &PixelSelectorAll::normalization() const noexcept {
  return _selectors.front().normalization();
}
inline MatrixUnit PixelSelectorAll::unit() const noexcept { return _selectors.front().unit(); }
inline std::uint32_t PixelSelectorAll::resolution() const noexcept {
  return _selectors.front().resolution();
}

inline const BinTable &PixelSelectorAll::bins() const noexcept {
  assert(_bins);
  return *_bins;
}
inline std::shared_ptr<const BinTable> PixelSelectorAll::bins_ptr() const noexcept { return _bins; }

inline const balancing::Weights &PixelSelectorAll::weights() const {
  if (!_weight_cache) {
    throw std::runtime_error(
        "PixelSelectorAll::weights() was called on an instance of PixelSelectorAll with null "
        "_weight_cache");
  }
  auto weights = _weight_cache->get_or_init(0, normalization());
  if (!weights->empty()) {
    return *weights;
  }

  if (normalization() == balancing::Method::NONE()) {
    *weights = balancing::Weights{1.0, bins().size(), balancing::Weights::Type::DIVISIVE};
    return *weights;
  }

  std::vector<double> buff(bins().size(), std::numeric_limits<double>::quiet_NaN());
  std::for_each(_selectors.begin(), _selectors.end(), [&](const PixelSelector &sel) {
    if (sel.is_intra()) {
      const auto &chrom_weights = sel.weights1();
      const auto offset = static_cast<std::ptrdiff_t>(bins().at(sel.chrom1()).id());
      std::copy(chrom_weights.begin(balancing::Weights::Type::DIVISIVE),
                chrom_weights.end(balancing::Weights::Type::DIVISIVE), buff.begin() + offset);
    }
  });

  *weights = balancing::Weights{std::move(buff), balancing::Weights::Type::DIVISIVE};
  return *weights;
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::Pair::operator<(const Pair &other) const noexcept {
  return first < other.first;
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::Pair::operator>(const Pair &other) const noexcept {
  return first > other.first;
}

template <typename N>
inline PixelSelectorAll::iterator<N>::iterator(const PixelSelectorAll &selector, bool sorted)
    : _selectors(std::make_shared<SelectorQueue>()),
      _active_selectors(std::make_shared<SelectorQueue>()),
      _its(std::make_shared<ItPQueue>()),
      _sorted(sorted),
      _buff(std::make_shared<std::vector<ThinPixel<N>>>()) {
  if (selector._selectors.empty()) {
    *this = iterator<N>{};
    return;
  }
  std::for_each(selector._selectors.begin(), selector._selectors.end(),
                [&](const PixelSelector &sel) { _selectors->push(&sel); });

  _chrom1_id = _selectors->front()->chrom1().id();
  init_iterators();
  read_next_chunk();
}

template <typename N>
inline PixelSelectorAll::iterator<N>::iterator(const iterator &other)
    : _selectors(other._selectors ? std::make_shared<SelectorQueue>(*other._selectors) : nullptr),
      _active_selectors(other._active_selectors
                            ? std::make_shared<SelectorQueue>(*other._active_selectors)
                            : nullptr),
      _its(other._its ? std::make_shared<ItPQueue>(*other._its) : nullptr),
      _sorted(other._sorted),
      _chrom1_id(other._chrom1_id),
      _buff(other._buff ? std::make_shared<std::vector<ThinPixel<N>>>(*other._buff) : nullptr),
      _i(other._i) {}

template <typename N>
inline auto PixelSelectorAll::iterator<N>::operator=(const iterator<N> &other) -> iterator & {
  if (this == &other) {
    return *this;
  }

  _selectors = other._selectors ? std::make_shared<SelectorQueue>(*other._selectors) : nullptr;
  _active_selectors =
      other._active_selectors ? std::make_shared<SelectorQueue>(*other._active_selectors) : nullptr;
  _its = other._its ? std::make_shared<ItPQueue>(*other._its) : nullptr;
  _sorted = other._sorted;
  _chrom1_id = other._chrom1_id;
  _buff = other._buff ? std::make_shared<std::vector<ThinPixel<N>>>(*other._buff) : nullptr;
  _i = other._i;

  return *this;
}

template <typename N>
inline bool PixelSelectorAll::iterator<N>::operator==(const iterator<N> &other) const noexcept {
  if (!_buff || !other._buff) {
    return _buff == other._buff;
  }

  assert(_i < _buff->size());
  assert(other._i < other._buff->size());
  return (*_buff)[_i] == (*other._buff)[other._i];
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

  if (_selectors.use_count() != 1) {
    _selectors = std::make_shared<SelectorQueue>(*_selectors);
  }

  if (_active_selectors.use_count() != 1) {
    _active_selectors = std::make_shared<SelectorQueue>(*_active_selectors);
  }

  while (!_active_selectors->empty()) {
    _active_selectors->front()->clear_cache();
    _active_selectors->pop();
  }

  while (!_selectors->empty() && _selectors->front()->chrom1().id() == _chrom1_id) {
    const auto *sel = _selectors->front();
    _selectors->pop();
    auto first = sel->begin<N>(_sorted);
    auto last = sel->end<N>();
    if (first != last) {
      _its->emplace(Pair{first, last});
      _active_selectors->emplace(sel);
    }
  }
}

template <typename N>
inline void PixelSelectorAll::iterator<N>::read_next_chunk() {
  constexpr std::size_t chunk_size = 100'000;
  if (_selectors->empty() && _its->empty()) {
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
    _buff = std::make_shared<std::vector<ThinPixel<N>>>();
  }
  _buff->clear();
  _i = 0;

  if (!_sorted) {
    while (first != last && _buff->size() != chunk_size) {
      _buff->push_back(*first);
      ++first;
    }
  } else {
    const auto bin1_id = first->bin1_id;
    while (first != last && first->bin1_id == bin1_id) {
      _buff->push_back(*first);
      ++first;
    }
  }

  _its->emplace(Pair{first, last});
}

}  // namespace hictk::hic
