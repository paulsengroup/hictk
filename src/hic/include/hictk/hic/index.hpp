// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>

#include "hictk/chromosome.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/hic_file_stream.hpp"

namespace hictk::hic::internal {

struct BlockIndex {
  std::size_t id{null_id};              // NOLINT
  std::size_t file_offset{};            // NOLINT
  std::size_t compressed_size_bytes{};  // NOLINT

  std::size_t first_row{};
  std::size_t last_row{};
  std::size_t first_col{};
  std::size_t last_col{};

  static constexpr auto null_id = (std::numeric_limits<std::size_t>::max)();

  constexpr explicit operator bool() const noexcept;
  friend constexpr bool operator<(const BlockIndex& a, const BlockIndex& b) noexcept;
  friend constexpr bool operator==(const BlockIndex& a, const BlockIndex& b) noexcept;
  friend constexpr bool operator!=(const BlockIndex& a, const BlockIndex& b) noexcept;

  friend constexpr bool operator<(const BlockIndex& a, std::size_t b_id) noexcept;
  friend constexpr bool operator==(const BlockIndex& a, std::size_t b_id) noexcept;
  friend constexpr bool operator!=(const BlockIndex& a, std::size_t b_id) noexcept;

  friend constexpr bool operator<(std::size_t a_id, const BlockIndex& b) noexcept;
  friend constexpr bool operator==(std::size_t a_id, const BlockIndex& b) noexcept;
  friend constexpr bool operator!=(std::size_t a_id, const BlockIndex& b) noexcept;
};

struct BlockIndexCmp {
  using is_transparent = void;

  constexpr bool operator()(const BlockIndex& a, const BlockIndex& b) const noexcept;
  constexpr bool operator()(const BlockIndex& a, std::size_t b_id) const noexcept;
  constexpr bool operator()(std::size_t a_id, const BlockIndex& b) const noexcept;
};

// Map coordinates (bp) to block IDs
class Index {
  // map block_ids to file offsets
  const phmap::btree_set<BlockIndex, BlockIndexCmp> _block_map{};
  std::int32_t _version{};
  std::size_t _block_bin_count{};
  std::size_t _block_column_count{};  // columns of blocks per matrix?
  double _sum_count{};                // sum

  Chromosome _chrom1{};
  Chromosome _chrom2{};

 public:
  static constexpr auto npos = (std::numeric_limits<std::size_t>::max)();

  Index() = default;
  Index(Chromosome chrom1_, Chromosome chrom2_, phmap::btree_set<BlockIndex, BlockIndexCmp> blocks_,
        std::int32_t version_, std::size_t block_bin_count_, std::size_t block_column_count_,
        double sum_count_);

  [[nodiscard]] const Chromosome& chrom1() const noexcept;
  [[nodiscard]] const Chromosome& chrom2() const noexcept;
  [[nodiscard]] bool is_intra() const noexcept;
  [[nodiscard]] constexpr double matrix_sum() const noexcept;
  [[nodiscard]] constexpr std::size_t block_bin_count() const noexcept;
  [[nodiscard]] constexpr std::size_t block_column_count() const noexcept;

  std::vector<BlockIndex> map_2d_query_to_blocks(const PixelCoordinates& coords1,
                                                 const PixelCoordinates& coords2);
  void map_2d_query_to_blocks(const PixelCoordinates& coords1, const PixelCoordinates& coords2,
                              std::vector<BlockIndex>& buffer);

  const BlockIndex& at(std::size_t id) const;

 private:
  void _map_2d_query_to_blocks(const PixelCoordinates& coords1, const PixelCoordinates& coords2,
                               std::vector<BlockIndex>& buffer);
  void _map_2d_query_to_blocks_intra_v9plus(const PixelCoordinates& coords1,
                                            const PixelCoordinates& coords2,
                                            std::vector<BlockIndex>& buffer);
};

}  // namespace hictk::hic::internal

#include "../../../index_impl.hpp"
