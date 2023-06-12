// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <memory>
#include <type_traits>
#include <utility>

namespace hictk::hic::internal {

constexpr BlockIndex::operator bool() const noexcept {
  return id != null_id && compressed_size_bytes != 0;
}

constexpr bool operator<(const BlockIndex &a, const BlockIndex &b) noexcept { return a < b.id; }
constexpr bool operator==(const BlockIndex &a, const BlockIndex &b) noexcept { return a == b.id; }
constexpr bool operator!=(const BlockIndex &a, const BlockIndex &b) noexcept { return !(a == b); }

constexpr bool operator<(const BlockIndex &a, std::size_t b_id) noexcept { return a.id < b_id; }
constexpr bool operator==(const BlockIndex &a, std::size_t b_id) noexcept { return a.id == b_id; }
constexpr bool operator!=(const BlockIndex &a, std::size_t b_id) noexcept { return !(a == b_id); }

constexpr bool operator<(std::size_t a_id, const BlockIndex &b) noexcept { return a_id < b.id; }
constexpr bool operator==(std::size_t a_id, const BlockIndex &b) noexcept { return a_id == b.id; }
constexpr bool operator!=(std::size_t a_id, const BlockIndex &b) noexcept { return !(a_id == b); }

constexpr bool BlockIndexCmp::operator()(const BlockIndex &a, const BlockIndex &b) const noexcept {
  return a < b;
}
constexpr bool BlockIndexCmp::operator()(const BlockIndex &a, std::size_t b_id) const noexcept {
  return a < b_id;
}
constexpr bool BlockIndexCmp::operator()(std::size_t a_id, const BlockIndex &b) const noexcept {
  return a_id < b;
}

inline Index::Index(Chromosome chrom1_, Chromosome chrom2_,
                    phmap::btree_set<BlockIndex, BlockIndexCmp> blocks_, std::int32_t version_,
                    std::size_t block_bin_count_, std::size_t block_column_count_,
                    double sum_count_)
    : _block_map(std::move(blocks_)),
      _version(version_),
      _block_bin_count(block_bin_count_),
      _block_column_count(block_column_count_),
      _sum_count(sum_count_),
      _chrom1(std::move(chrom1_)),
      _chrom2(std::move(chrom2_)) {
  if (_block_bin_count == 0) {
    throw std::runtime_error("index is corrupted: blockBinCount=0.");
  }
  if (_block_column_count == 0) {
    throw std::runtime_error("index is corrupted: blockColumnCount=0.");
  }
}

inline const Chromosome &Index::chrom1() const noexcept { return _chrom1; }
inline const Chromosome &Index::chrom2() const noexcept { return _chrom2; }
inline bool Index::is_intra() const noexcept { return _chrom1 == _chrom2; }

constexpr double Index::matrix_sum() const noexcept { return _sum_count; }

constexpr std::size_t Index::block_bin_count() const noexcept { return _block_bin_count; }

constexpr std::size_t Index::block_column_count() const noexcept { return _block_column_count; }

inline std::vector<BlockIndex> Index::map_2d_query_to_blocks(const PixelCoordinates &coords1,
                                                             const PixelCoordinates &coords2) {
  std::vector<BlockIndex> buffer{};
  map_2d_query_to_blocks(coords1, coords2, buffer);
  return buffer;
}

inline const BlockIndex &Index::at(std::size_t id) const {
  auto match = _block_map.find(id);
  if (match == _block_map.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("unable to find block #{}: out of range"), id));
  }
  return *match;
}

inline void Index::map_2d_query_to_blocks(const PixelCoordinates &coords1,
                                          const PixelCoordinates &coords2,
                                          std::vector<BlockIndex> &buffer) {
  assert(coords1.is_intra());
  assert(coords2.is_intra());

  const auto is_intra = coords1.bin1.chrom() == coords2.bin1.chrom();
  if (_version < 9 || is_intra) {
    return _map_2d_query_to_blocks(coords1, coords2, buffer);
  }
  return _map_2d_query_to_blocks_intra_v9plus(coords1, coords2, buffer);
}

inline void Index::_map_2d_query_to_blocks(const hictk::PixelCoordinates &coords1,
                                           const hictk::PixelCoordinates &coords2,
                                           std::vector<BlockIndex> &buffer) {
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
  phmap::btree_set<BlockIndex> tmp_buffer{};
  // first check the upper triangular matrix_type
  for (auto row = row1; row <= row2; ++row) {
    for (auto col = col1; col <= col2; ++col) {
      const auto id1 = row * _block_column_count + col;
      auto match = _block_map.find(id1);
      if (match != _block_map.end()) {
        auto block = *match;
        block.first_row = bin1;
        block.last_row = bin2;
        block.first_col = bin3;
        block.last_col = bin4;

        tmp_buffer.emplace(block);
      }

      if (checkLowerLeftTri) {
        const auto id2 = col * _block_column_count + row;
        match = _block_map.find(id2);
        if (match != _block_map.end()) {
          auto block = *match;
          block.first_row = bin1;
          block.last_row = bin2;
          block.first_col = bin3;
          block.last_col = bin4;
          tmp_buffer.emplace(block);
        }
      }
    }
  }

  buffer.resize(tmp_buffer.size());
  std::copy(tmp_buffer.begin(), tmp_buffer.end(), buffer.begin());
}

inline void Index::_map_2d_query_to_blocks_intra_v9plus(const hictk::PixelCoordinates &coords1,
                                                        const hictk::PixelCoordinates &coords2,
                                                        std::vector<BlockIndex> &buffer) {
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
  phmap::btree_set<BlockIndex> block_ids{};
  for (auto depth = nearerDepth; depth <= furtherDepth; ++depth) {
    for (auto pad = translatedLowerPAD; pad <= translatedHigherPAD; ++pad) {
      const auto id = depth * _block_column_count + pad;
      auto match = _block_map.find(id);
      if (match != _block_map.end()) {
        auto block = *match;
        block.first_row = bin1;
        block.last_row = bin2;
        block.first_col = bin3;
        block.last_col = bin3;
        block_ids.emplace(block);
      }
    }
  }
  buffer.resize(block_ids.size());
  std::copy(block_ids.begin(), block_ids.end(), buffer.begin());
}

}  // namespace hictk::hic::internal
