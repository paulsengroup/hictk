// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>

namespace hictk::hic::internal {

constexpr bool BlockIndex::GridCoordinates::operator==(
    const GridCoordinates &other) const noexcept {
  return i1 == other.i1 && i2 == other.i2;
}

constexpr bool BlockIndex::GridCoordinates::operator!=(
    const GridCoordinates &other) const noexcept {
  return !(*this == other);
}

constexpr bool BlockIndex::GridCoordinates::operator<(const GridCoordinates &other) const noexcept {
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

constexpr double Index::matrix_sum() const noexcept { return _sum_count; }

constexpr std::size_t Index::block_bin_count() const noexcept { return _block_bin_count; }

constexpr std::size_t Index::block_column_count() const noexcept { return _block_column_count; }

}  // namespace hictk::hic::internal
