// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>

#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

class BlockIndex {
 public:
  struct GridCoordinates {
    std::size_t row;
    std::size_t col;

    constexpr bool operator==(const GridCoordinates& other) const noexcept;
    constexpr bool operator!=(const GridCoordinates& other) const noexcept;
    constexpr bool operator<(const GridCoordinates& other) const noexcept;
  };

  std::size_t _id{null_id};              // NOLINT
  std::size_t _file_offset{};            // NOLINT
  std::size_t _compressed_size_bytes{};  // NOLINT
  GridCoordinates _coords{};             // NOLINT

  static constexpr auto null_id = (std::numeric_limits<std::size_t>::max)();

  constexpr BlockIndex() = default;
  constexpr BlockIndex(std::size_t id_, std::size_t file_offset_,
                       std::size_t compressed_size_bytes_, std::size_t block_column_count) noexcept;

  [[nodiscard]] constexpr std::size_t id() const noexcept;
  [[nodiscard]] constexpr std::size_t file_offset() const noexcept;
  [[nodiscard]] constexpr std::size_t compressed_size_bytes() const noexcept;
  [[nodiscard]] constexpr auto coords() const noexcept -> const GridCoordinates&;

  constexpr explicit operator bool() const noexcept;
  constexpr bool operator==(const BlockIndex& other) const noexcept;
  constexpr bool operator!=(const BlockIndex& other) const noexcept;
  constexpr bool operator<(const BlockIndex& other) const noexcept;
  constexpr bool operator==(const BlockIndex::GridCoordinates& coords_) const noexcept;
  constexpr bool operator!=(const BlockIndex::GridCoordinates& coords_) const noexcept;
  constexpr bool operator<(const BlockIndex::GridCoordinates& coords_) const noexcept;
};

struct BlockIndexCmp {
  using is_transparent = void;

  constexpr bool operator()(const BlockIndex& a, const BlockIndex& b) const noexcept;
  constexpr bool operator()(const BlockIndex& a,
                            const BlockIndex::GridCoordinates& b_coords) const noexcept;
  constexpr bool operator()(const BlockIndex::GridCoordinates& a_coords,
                            const BlockIndex& b) const noexcept;
};

// Map coordinates (bp) to block IDs
class Index {
  using BlockIndexMap = phmap::btree_set<BlockIndex, BlockIndexCmp>;
  // map block_ids to file offsets
  const BlockIndexMap _block_map{};
  std::int32_t _version{};
  std::size_t _block_bin_count{};
  std::size_t _block_column_count{};  // columns of blocks per matrix?
  double _sum_count{};                // sum

  MatrixUnit _unit{};
  std::uint32_t _resolution{};
  Chromosome _chrom1{};
  Chromosome _chrom2{};

 public:
  static constexpr auto npos = (std::numeric_limits<std::size_t>::max)();

  Index() = default;
  Index(Chromosome chrom1_, Chromosome chrom2_, MatrixUnit unit_, std::uint32_t resolution_,
        std::int32_t version_, std::size_t block_bin_count_, std::size_t block_column_count_,
        double sum_count_, BlockIndexMap blocks_);

  [[nodiscard]] MatrixUnit unit() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] const Chromosome& chrom1() const noexcept;
  [[nodiscard]] const Chromosome& chrom2() const noexcept;
  [[nodiscard]] bool is_intra() const noexcept;
  [[nodiscard]] constexpr double matrix_sum() const noexcept;
  [[nodiscard]] constexpr std::size_t block_bin_count() const noexcept;
  [[nodiscard]] constexpr std::size_t block_column_count() const noexcept;

  [[nodiscard]] auto begin() const noexcept -> BlockIndexMap::const_iterator;
  [[nodiscard]] auto end() const noexcept -> BlockIndexMap::const_iterator;
  [[nodiscard]] auto cbegin() const noexcept -> BlockIndexMap::const_iterator;
  [[nodiscard]] auto cend() const noexcept -> BlockIndexMap::const_iterator;

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;

  [[nodiscard]] Index subset(const PixelCoordinates& coords1,
                             const PixelCoordinates& coords2) const;

  [[nodiscard]] auto map_2d_query_to_blocks(const PixelCoordinates& coords1,
                                            const PixelCoordinates& coords2) const -> BlockIndexMap;

  [[nodiscard]] const BlockIndex& at(std::size_t row, std::size_t col) const;

 private:
  void _map_2d_query_to_blocks(const PixelCoordinates& coords1, const PixelCoordinates& coords2,
                               BlockIndexMap& buffer) const;
  void _map_2d_query_to_blocks_intra_v9plus(const PixelCoordinates& coords1,
                                            const PixelCoordinates& coords2,
                                            BlockIndexMap& buffer) const;
};

}  // namespace hictk::hic::internal

#include "../../../index_impl.hpp"
