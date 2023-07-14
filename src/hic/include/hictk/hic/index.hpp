// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

class BlockIndex {
 public:
  struct GridCoordinates {
    std::size_t row;  // NOLINT
    std::size_t col;  // NOLINT

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
  constexpr bool operator==(std::size_t id_) const noexcept;
  constexpr bool operator!=(std::size_t id_) const noexcept;
};

// Map coordinates (bp) to block IDs
class Index {
 public:
  using BlkIdxBuffer = std::vector<BlockIndex>;
  using iterator = BlkIdxBuffer::const_iterator;
  using const_iterator = BlkIdxBuffer::const_iterator;

 private:
  // map block_ids to file offsets
  BlkIdxBuffer _buffer{};
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
        double sum_count_, BlkIdxBuffer blocks_);

  [[nodiscard]] MatrixUnit unit() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] const Chromosome& chrom1() const noexcept;
  [[nodiscard]] const Chromosome& chrom2() const noexcept;
  [[nodiscard]] bool is_intra() const noexcept;
  [[nodiscard]] constexpr double matrix_sum() const noexcept;
  [[nodiscard]] constexpr std::size_t block_bin_count() const noexcept;
  [[nodiscard]] constexpr std::size_t block_column_count() const noexcept;

  [[nodiscard]] auto begin() const noexcept -> const_iterator;
  [[nodiscard]] auto end() const noexcept -> const_iterator;
  [[nodiscard]] auto cbegin() const noexcept -> const_iterator;
  [[nodiscard]] auto cend() const noexcept -> const_iterator;

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;

  [[nodiscard]] auto find_overlaps(const Bin& bin1, const PixelCoordinates& coords2) const
      -> std::pair<const_iterator, const_iterator>;

  [[nodiscard]] const BlockIndex& at(std::size_t row, std::size_t col) const;

 private:
  [[nodiscard]] auto generate_block_list(std::size_t bin1, std::size_t bin2, std::size_t bin3,
                                         std::size_t bin4) const
      -> std::pair<const_iterator, const_iterator>;
  [[nodiscard]] auto generate_block_list_intra_v9plus(std::size_t bin1, std::size_t bin2,
                                                      std::size_t bin3, std::size_t bin4) const
      -> std::pair<const_iterator, const_iterator>;
};

}  // namespace hictk::hic::internal

template <>
struct std::hash<hictk::hic::internal::BlockIndex> {
  inline std::size_t operator()(hictk::hic::internal::BlockIndex const& b) const noexcept {
    return std::hash<std::size_t>{}(b.id());
  }
};
#include "../../../index_impl.hpp"
