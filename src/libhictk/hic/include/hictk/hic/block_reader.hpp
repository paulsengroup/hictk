// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/binary_buffer.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

class HiCBlockReader {
  std::shared_ptr<HiCFileReader> _hfs{};
  std::shared_ptr<BlockCache> _blk_cache{};
  // We need the entire bin table in order to map pixels to abs bin ids
  std::shared_ptr<const BinTable> _bins{};
  Index _index{};

  BinaryBuffer _bbuffer{};
  std::vector<ThinPixel<float>> _tmp_buffer{};

 public:
  HiCBlockReader() = default;
  HiCBlockReader(std::shared_ptr<HiCFileReader> hfs, Index master_index,
                 std::shared_ptr<const BinTable> bins_,
                 std::shared_ptr<BlockCache> block_cache_) noexcept;

  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] const Chromosome& chrom1() const noexcept;
  [[nodiscard]] const Chromosome& chrom2() const noexcept;
  [[nodiscard]] const BinTable& bins() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;
  [[nodiscard]] const Index& index() const noexcept;

  [[nodiscard]] double sum() const noexcept;
  [[nodiscard]] double avg() const;

  [[nodiscard]] std::shared_ptr<const InteractionBlock> read_v6(const Chromosome& chrom1,
                                                                const Chromosome& chrom2,
                                                                const BlockIndex& idx,
                                                                bool cache_block = true);
  [[nodiscard]] std::shared_ptr<const InteractionBlock> read(const Chromosome& chrom1,
                                                             const Chromosome& chrom2,
                                                             const BlockIndex& idx,
                                                             bool cache_block = true);
  [[nodiscard]] std::size_t read_size(const Chromosome& chrom1, const Chromosome& chrom2,
                                      const BlockIndex& idx);
  void evict(const InteractionBlock& blk);
  void evict(const Chromosome& chrom1, const Chromosome& chrom2, const BlockIndex& idx);
  void clear() noexcept;

  [[nodiscard]] std::size_t cache_size() const noexcept;

 private:
  static void read_dispatcher_type1_block(bool i16Bin1, bool i16Bin2, bool i16Counts,
                                          std::int32_t bin1Offset, std::int32_t bin2Offset,
                                          BinaryBuffer& src,
                                          std::vector<ThinPixel<float>>& dest) noexcept;
  template <typename Bin1Type, typename Bin2Type, typename CountType>
  static void read_type1_block(std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer& src,
                               std::vector<ThinPixel<float>>& dest) noexcept;

  template <typename CountType>
  static void read_type2_block(std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer& src,
                               std::vector<ThinPixel<float>>& dest) noexcept;
};

}  // namespace hictk::hic::internal

#include "./impl/block_reader_impl.hpp"  // NOLINT
