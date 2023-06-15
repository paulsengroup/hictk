// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/numeric_utils.hpp"

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

constexpr bool BlockIndex::operator==(std::size_t id_) const noexcept { return _id == id_; }

constexpr bool BlockIndex::operator!=(std::size_t id_) const noexcept { return !(*this == id_); }

inline std::size_t BlockIndexHasher::operator()(const BlockIndex &b) const noexcept {
  return (*this)(b.id());
}
inline std::size_t BlockIndexHasher::operator()(std::size_t id) const noexcept {
  return std::hash<std::size_t>{}(id);
}

constexpr bool BlockIndexEq::operator()(const BlockIndex &a, const BlockIndex &b) const noexcept {
  return a == b;
}
constexpr bool BlockIndexEq::operator()(std::size_t a_id, const BlockIndex &b) const noexcept {
  return a_id == b.id();
}
constexpr bool BlockIndexEq::operator()(const BlockIndex &a, std::size_t b_id) const noexcept {
  return a.id() == b_id;
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

inline std::vector<BlockIndex> Index::find_overlaps(const PixelCoordinates &coords1,
                                                    const PixelCoordinates &coords2) const {
  assert(coords1.is_intra());
  assert(coords2.is_intra());

  std::vector<BlockIndex> buffer{};
  map_2d_query_to_blocks(coords1, coords2, buffer);
  return buffer;
}

inline const BlockIndex &Index::at(std::size_t row, std::size_t col) const {
  const auto block_id = (col * block_column_count()) + row;
  auto match = _block_map.find(block_id);
  if (match == _block_map.end()) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("unable to find block {}{}: out of range"), row, col));
  }
  return *match;
}

inline void Index::generate_block_list(std::size_t bin1, std::size_t bin2, std::size_t bin3,
                                       std::size_t bin4, bool is_intra) const {
  const auto col1 = bin1 / _block_bin_count;
  const auto col2 = (bin2 + 1) / _block_bin_count;
  const auto row1 = bin3 / _block_bin_count;
  const auto row2 = (bin4 + 1) / _block_bin_count;

  // check region part that overlaps with lower left triangle but only if intrachromosomal
  for (auto row = row1; row <= row2; ++row) {
    for (auto col = col1; col <= col2; ++col) {
      const auto block_id = (row * block_column_count()) + col;
      const auto match = _block_map.find(block_id);
      if (match != _block_map.end()) {
        _tmp_buffer.emplace(*match);
      }
    }
  }

  if (is_intra) {
    std::swap(bin1, bin3);
    std::swap(bin3, bin4);
    generate_block_list(bin1, bin2, bin3, bin4, false);
  }
}

inline void Index::generate_block_list_intra_v9plus(std::size_t bin1, std::size_t bin2,
                                                    std::size_t bin3, std::size_t bin4) const {
  const auto translatedLowerPAD = (bin1 + bin3) / 2 / _block_bin_count;
  const auto translatedHigherPAD = (bin2 + bin4) / 2 / _block_bin_count + 1;
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
    return (std::min)(translatedNearerDepth, translatedFurtherDepth);
  }();

  // +1; integer divide rounds down
  const auto furtherDepth = (std::max)(translatedNearerDepth, translatedFurtherDepth) + 1;
  for (auto depth = nearerDepth; depth <= furtherDepth; ++depth) {
    for (auto pad = translatedLowerPAD; pad <= translatedHigherPAD; ++pad) {
      const auto block_id = (depth * block_column_count()) + pad;
      auto match = _block_map.find(block_id);
      if (match != _block_map.end()) {
        _tmp_buffer.emplace(*match);
      }
    }
  }
}

inline void Index::map_2d_query_to_blocks(const hictk::PixelCoordinates &coords1,
                                          const hictk::PixelCoordinates &coords2,
                                          std::vector<BlockIndex> &buffer) const {
  assert(coords1.bin1.chrom() == _chrom1 || coords1.bin1.chrom() == _chrom2);
  assert(coords2.bin1.chrom() == _chrom1 || coords2.bin1.chrom() == _chrom2);

  auto bin1 = coords1.bin1.rel_id();
  auto bin2 = coords1.bin2.rel_id() + 1;
  auto bin3 = coords2.bin1.rel_id();
  auto bin4 = coords2.bin2.rel_id() + 1;

  const auto is_intra = coords1.bin1.chrom() == coords2.bin1.chrom();

  _tmp_buffer.clear();
  if (_version > 8 && is_intra) {
    generate_block_list_intra_v9plus(bin1, bin2, bin3, bin4);
  } else {
    generate_block_list(bin1, bin2, bin3, bin4, is_intra);
  }

  buffer.resize(_tmp_buffer.size());
  std::move(_tmp_buffer.begin(), _tmp_buffer.end(), buffer.begin());
  std::sort(buffer.begin(), buffer.end());
}

}  // namespace hictk::hic::internal
