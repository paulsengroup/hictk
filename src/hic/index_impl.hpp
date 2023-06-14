// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace hictk::hic::internal {

constexpr bool BlockIndex::GridCoordinates::operator==(
    const BlockIndex::GridCoordinates &other) const noexcept {
  return row == other.row && col == other.col;
}

constexpr bool BlockIndex::GridCoordinates::operator!=(
    const BlockIndex::GridCoordinates &other) const noexcept {
  return !(*this == other);
}

constexpr bool BlockIndex::GridCoordinates::operator<(
    const BlockIndex::GridCoordinates &other) const noexcept {
  if (row == other.row) {
    return col < other.col;
  }
  return row < other.row;
}

constexpr BlockIndex::BlockIndex(std::size_t id_, std::size_t file_offset_,
                                 std::size_t compressed_size_bytes_,
                                 std::size_t block_column_count) noexcept
    : _id(id_),
      _file_offset(file_offset_),
      _compressed_size_bytes(compressed_size_bytes_),
      _coords({_id % block_column_count, _id / block_column_count}) {}

constexpr std::size_t BlockIndex::id() const noexcept { return _id; }
constexpr std::size_t BlockIndex::file_offset() const noexcept { return _file_offset; }
constexpr std::size_t BlockIndex::compressed_size_bytes() const noexcept {
  return _compressed_size_bytes;
}
constexpr auto BlockIndex::coords() const noexcept -> const GridCoordinates & { return _coords; }

constexpr BlockIndex::operator bool() const noexcept {
  return _id != null_id && _compressed_size_bytes != 0;
}

constexpr bool BlockIndex::operator==(const BlockIndex &other) const noexcept {
  return _id == other._id;
}

constexpr bool BlockIndex::operator!=(const BlockIndex &other) const noexcept {
  return !(*this == other);
}

constexpr bool BlockIndex::operator<(const BlockIndex &other) const noexcept {
  return _coords < other._coords;
}

constexpr bool BlockIndex::operator==(const BlockIndex::GridCoordinates &coords_) const noexcept {
  return _coords == coords_;
}

constexpr bool BlockIndex::operator!=(const BlockIndex::GridCoordinates &coords_) const noexcept {
  return !(*this == coords_);
}

constexpr bool BlockIndex::operator<(const BlockIndex::GridCoordinates &coords_) const noexcept {
  return _coords < coords_;
}

constexpr bool BlockIndexCmp::operator()(const BlockIndex &a, const BlockIndex &b) const noexcept {
  return a.coords() < b.coords();
}

constexpr bool BlockIndexCmp::operator()(
    const BlockIndex &a, const BlockIndex::GridCoordinates &b_coords) const noexcept {
  return a < b_coords;
}
constexpr bool BlockIndexCmp::operator()(const BlockIndex::GridCoordinates &a_coords,
                                         const BlockIndex &b) const noexcept {
  return a_coords < b._coords;
}

inline Index::Index(Chromosome chrom1_, Chromosome chrom2_, MatrixUnit unit_,
                    std::uint32_t resolution_, std::int32_t version_, std::size_t block_bin_count_,
                    std::size_t block_column_count_, double sum_count_, BlockIndexMap blocks_)
    : _block_map(std::move(blocks_)),
      _version(version_),
      _block_bin_count(block_bin_count_),
      _block_column_count(block_column_count_),
      _sum_count(sum_count_),
      _unit(unit_),
      _resolution(resolution_),
      _chrom1(std::move(chrom1_)),
      _chrom2(std::move(chrom2_)) {}

inline MatrixUnit Index::unit() const noexcept { return _unit; }
inline std::uint32_t Index::resolution() const noexcept { return _resolution; }
inline const Chromosome &Index::chrom1() const noexcept { return _chrom1; }
inline const Chromosome &Index::chrom2() const noexcept { return _chrom2; }
inline bool Index::is_intra() const noexcept { return _chrom1 == _chrom2; }

constexpr double Index::matrix_sum() const noexcept { return _sum_count; }

constexpr std::size_t Index::block_bin_count() const noexcept { return _block_bin_count; }

constexpr std::size_t Index::block_column_count() const noexcept { return _block_column_count; }

inline auto Index::begin() const noexcept -> BlockIndexMap::const_iterator {
  return _block_map.begin();
}
inline auto Index::end() const noexcept -> BlockIndexMap::const_iterator {
  return _block_map.end();
}
inline auto Index::cbegin() const noexcept -> BlockIndexMap::const_iterator {
  return _block_map.cbegin();
}
inline auto Index::cend() const noexcept -> BlockIndexMap::const_iterator {
  return _block_map.cend();
}

inline std::size_t Index::size() const noexcept { return _block_map.size(); }

inline bool Index::empty() const noexcept { return size() == 0; }  // NOLINT

inline Index Index::subset(const PixelCoordinates &coords1, const PixelCoordinates &coords2) const {
  return {
      chrom1(),
      chrom2(),
      unit(),
      resolution(),
      _version,
      block_bin_count(),
      block_column_count(),
      matrix_sum(),
      map_2d_query_to_blocks(coords1, coords2),
  };
}

inline const BlockIndex &Index::at(std::size_t row, std::size_t col) const {
  auto match = _block_map.find(BlockIndex::GridCoordinates{row, col});
  if (match == _block_map.end()) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("unable to find block {}{}: out of range"), row, col));
  }
  return *match;
}

inline auto Index::map_2d_query_to_blocks(const PixelCoordinates &coords1,
                                          const PixelCoordinates &coords2) const -> BlockIndexMap {
  assert(coords1.is_intra());
  assert(coords2.is_intra());

  BlockIndexMap buffer{};

  const auto is_intra = coords1.bin1.chrom() == coords2.bin1.chrom();
  if (_version < 9 || is_intra) {
    _map_2d_query_to_blocks(coords1, coords2, buffer);
  } else {
    _map_2d_query_to_blocks_intra_v9plus(coords1, coords2, buffer);
  }
  return buffer;
}

inline void Index::_map_2d_query_to_blocks(const hictk::PixelCoordinates &coords1,
                                           const hictk::PixelCoordinates &coords2,
                                           BlockIndexMap &buffer) const {
  assert(coords1.bin1.chrom() == _chrom1 || coords1.bin1.chrom() == _chrom2);
  assert(coords2.bin1.chrom() == _chrom1 || coords2.bin1.chrom() == _chrom2);

  auto bin1 = coords1.bin1.rel_id();
  auto bin2 = coords1.bin2.rel_id() + 1;
  auto bin3 = coords2.bin1.rel_id();
  auto bin4 = coords2.bin2.rel_id() + 1;

  const auto is_intra = coords1.bin1.chrom() == coords2.bin1.chrom();

  if (is_intra && bin1 > bin3) {
    std::swap(bin1, bin3);
    std::swap(bin2, bin4);
  }

  const auto col1 = bin1 / _block_bin_count;
  const auto col2 = (bin2 + 1) / _block_bin_count;
  const auto row1 = bin3 / _block_bin_count;
  const auto row2 = (bin4 + 1) / _block_bin_count;

  // check region part that overlaps with lower left triangle but only if intrachromosomal
  const auto checkLowerLeftTri = is_intra;
  buffer.clear();
  // first check the upper triangular matrix_type
  for (auto row = row1; row <= row2; ++row) {
    for (auto col = col1; col <= col2; ++col) {
      auto match = _block_map.find(BlockIndex::GridCoordinates{row, col});
      if (match != _block_map.end()) {
        buffer.emplace(*match);
      }

      if (checkLowerLeftTri) {
        match = _block_map.find(BlockIndex::GridCoordinates{col, row});
        if (match != _block_map.end()) {
          buffer.emplace(*match);
        }
      }
    }
  }
}

inline void Index::_map_2d_query_to_blocks_intra_v9plus(const hictk::PixelCoordinates &coords1,
                                                        const hictk::PixelCoordinates &coords2,
                                                        BlockIndexMap &buffer) const {
  // https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#grid-structure
  assert(coords1.bin1.chrom() == _chrom1 || coords1.bin1.chrom() == _chrom2);
  assert(coords2.bin1.chrom() == _chrom1 || coords2.bin1.chrom() == _chrom2);

  auto bin1 = coords1.bin1.rel_id();
  auto bin2 = coords1.bin2.rel_id() + 1;
  auto bin3 = coords2.bin1.rel_id();
  auto bin4 = coords2.bin2.rel_id() + 1;

  const auto is_intra = coords1.bin1.chrom() == coords2.bin1.chrom();

  if (is_intra && bin1 > bin3) {
    std::swap(bin1, bin3);
    std::swap(bin2, bin4);
  }

  const auto translatedLowerPAD = (bin1 + bin3) / 2 / _block_bin_count;
  const auto translatedHigherPAD = (bin2 + bin4) / 2 / _block_column_count + 1;
  const auto translatedNearerDepth =
      static_cast<std::size_t>(std::log2(1.0 + double(hictk::internal::abs_diff(bin1, bin4)) /
                                                   std::sqrt(2.0) / double(_block_bin_count)));
  const auto translatedFurtherDepth =
      static_cast<std::size_t>(std::log2(1.0 + double(hictk::internal::abs_diff(bin2, bin3)) /
                                                   std::sqrt(2.0) / double(_block_bin_count)));

  // code above assumes diagonal; but we could be below diagonal
  const auto nearerDepth = [&]() -> std::size_t {
    if ((bin1 > bin4 && bin2 < bin3) || (bin2 > bin3 && bin1 < bin4)) {
      return 0;
    }
    return std::min(translatedNearerDepth, translatedFurtherDepth);
  }();

  // +1; integer divide rounds down
  const auto furtherDepth = std::max(translatedNearerDepth, translatedFurtherDepth) + 1;
  for (auto depth = nearerDepth; depth <= furtherDepth; ++depth) {
    for (auto pad = translatedLowerPAD; pad <= translatedHigherPAD; ++pad) {
      auto match = _block_map.find(BlockIndex::GridCoordinates{depth, pad});
      if (match != _block_map.end()) {
        auto block = *match;
        buffer.emplace(block);
      }
    }
  }
}

}  // namespace hictk::hic::internal
