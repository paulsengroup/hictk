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
  return i1 == other.i1 && i2 == other.i2;
}

constexpr bool BlockIndex::GridCoordinates::operator!=(
    const BlockIndex::GridCoordinates &other) const noexcept {
  return !(*this == other);
}

constexpr bool BlockIndex::GridCoordinates::operator<(
    const BlockIndex::GridCoordinates &other) const noexcept {
  if (i1 == other.i1) {
    return i2 < other.i2;
  }
  return i1 < other.i1;
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
  return _id < other._id;
}

constexpr bool BlockIndex::operator==(std::size_t id_) const noexcept { return _id == id_; }

constexpr bool BlockIndex::operator!=(std::size_t id_) const noexcept { return !(*this == id_); }

inline Index::Index(Chromosome chrom1_, Chromosome chrom2_, MatrixUnit unit_,
                    std::uint32_t resolution_, std::int32_t version_, std::size_t block_bin_count_,
                    std::size_t block_column_count_, double sum_count_, BlkIdxBuffer blocks_)
    : _buffer(std::move(blocks_)),
      _version(version_),
      _block_bin_count(block_bin_count_),
      _block_column_count(block_column_count_),
      _sum_count(sum_count_),
      _unit(unit_),
      _resolution(resolution_),
      _chrom1(std::move(chrom1_)),
      _chrom2(std::move(chrom2_)) {}

inline std::int32_t Index::version() const noexcept { return _version; }
inline MatrixUnit Index::unit() const noexcept { return _unit; }
inline std::uint32_t Index::resolution() const noexcept { return _resolution; }
inline const Chromosome &Index::chrom1() const noexcept { return _chrom1; }
inline const Chromosome &Index::chrom2() const noexcept { return _chrom2; }
inline bool Index::is_intra() const noexcept { return _chrom1 == _chrom2; }

constexpr double Index::matrix_sum() const noexcept { return _sum_count; }

constexpr std::size_t Index::block_bin_count() const noexcept { return _block_bin_count; }

constexpr std::size_t Index::block_column_count() const noexcept { return _block_column_count; }

inline auto Index::begin() const noexcept -> const_iterator { return _buffer.begin(); }
inline auto Index::end() const noexcept -> const_iterator { return _buffer.end(); }
inline auto Index::cbegin() const noexcept -> const_iterator { return _buffer.cbegin(); }
inline auto Index::cend() const noexcept -> const_iterator { return _buffer.cend(); }

inline std::size_t Index::size() const noexcept { return _buffer.size(); }

inline bool Index::empty() const noexcept { return size() == 0; }  // NOLINT

inline auto Index::find_overlaps(const PixelCoordinates &coords1,
                                 const PixelCoordinates &coords2) const -> Overlap {
  assert(coords1.is_intra());
  assert(coords2.is_intra());

  if (this->empty()) {
    return {};
  }

  assert(coords1.bin1.chrom() == _chrom1 || coords1.bin1.chrom() == _chrom2);
  assert(coords2.bin1.chrom() == _chrom1 || coords2.bin1.chrom() == _chrom2);

  auto bin1_id = coords1.bin1.rel_id();
  auto bin2_id = coords1.bin2.rel_id() + 1;
  auto bin3_id = coords2.bin1.rel_id();
  auto bin4_id = coords2.bin2.rel_id() + 1;

  const auto is_intra = coords1.bin1.chrom() == coords2.bin1.chrom();

  return _version > 8 && is_intra
             ? generate_block_list_intra_v9plus(bin1_id, bin2_id, bin3_id, bin4_id)
             : generate_block_list(bin1_id, bin2_id, bin3_id, bin4_id);
}

inline const BlockIndex &Index::at(std::size_t row, std::size_t col) const {
  const auto block_id = (col * block_column_count()) + row;
  const auto &[first, last] = std::equal_range(_buffer.begin(), _buffer.end(),
                                               BlockIndex{block_id, 0, 0, _block_column_count});
  if (first == last) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("unable to find block {}{}: out of range"), row, col));
  }
  assert(std::distance(first, last) == 1);
  return *first;
}

inline auto Index::generate_block_list(std::size_t bin1, std::size_t bin2, std::size_t bin3,
                                       std::size_t bin4) const -> Overlap {
  const auto col1 = bin1 / _block_bin_count;
  const auto col2 = (bin2 + 1) / _block_bin_count;
  const auto row1 = bin3 / _block_bin_count;
  const auto row2 = (bin4 + 1) / _block_bin_count;

  Overlap buffer{};
  for (auto row = row1; row <= row2; ++row) {
    for (auto col = col1; col <= col2; ++col) {
      const auto block_id = (row * block_column_count()) + col;
      const auto match = _buffer.find(BlockIndex{block_id, 0, 0, _block_column_count});
      if (match != _buffer.end()) {
        buffer.emplace_back(*match);
      }
    }
  }

  // Sort first by row, then by column
  std::sort(buffer.begin(), buffer.end(), [](const BlockIndex &b1, const BlockIndex &b2) {
    if (b1.coords().i1 != b2.coords().i1) {
      return b1.coords().i1 < b2.coords().i1;
    }
    return b1.coords().i2 < b2.coords().i2;
  });
  return buffer;
}

inline auto Index::generate_block_list_intra_v9plus(std::size_t bin1, std::size_t bin2,
                                                    std::size_t bin3, std::size_t bin4) const
    -> Overlap {
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

  const auto furtherDepth = (std::max)(translatedNearerDepth, translatedFurtherDepth) + 1;

  Overlap buffer{};
  for (auto pad = translatedLowerPAD; pad <= translatedHigherPAD; ++pad) {
    for (auto depth = nearerDepth; depth <= furtherDepth; ++depth) {
      const auto block_id = (depth * block_column_count()) + pad;
      const auto match = _buffer.find(BlockIndex{block_id, 0, 0, _block_column_count});
      if (match != _buffer.end()) {
        buffer.emplace_back(*match);
      }
    }
  }

  std::sort(buffer.begin(), buffer.end());
  buffer.erase(std::unique(buffer.begin(), buffer.end()), buffer.end());
  return buffer;
}

}  // namespace hictk::hic::internal
