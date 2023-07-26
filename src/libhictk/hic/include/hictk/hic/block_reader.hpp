// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>

#include "hictk/chromosome.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/index.hpp"

namespace hictk::hic::internal {

class BinaryBuffer {
  std::string _buffer{};
  std::size_t _i{};

 public:
  BinaryBuffer() = default;
  template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type* = nullptr>
  T read();

  // Return the offset of the underlying buffer. Useful for error checking
  [[nodiscard]] std::size_t operator()() const noexcept;

  // Reset and return ref to underlying buffer so that buff can be refilled
  std::string& reset() noexcept;
};

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
  HiCBlockReader(std::shared_ptr<HiCFileReader> hfs, const Index& master_index,
                 std::shared_ptr<const BinTable> bins_, std::shared_ptr<BlockCache> block_cache_);

  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] const Chromosome& chrom1() const noexcept;
  [[nodiscard]] const Chromosome& chrom2() const noexcept;
  [[nodiscard]] const BinTable& bins() const noexcept;
  [[nodiscard]] const Index& index() const noexcept;

  [[nodiscard]] double sum() const noexcept;
  [[nodiscard]] double avg() const;

  [[nodiscard]] std::shared_ptr<const InteractionBlock> read(const Chromosome& chrom1,
                                                             const Chromosome& chrom2,
                                                             const BlockIndex& idx);
  [[nodiscard]] std::size_t read_size(const Chromosome& chrom1, const Chromosome& chrom2,
                                      const BlockIndex& idx);
  void evict(const InteractionBlock& blk);
  void evict(const Chromosome& chrom1, const Chromosome& chrom2, const BlockIndex& idx);
  void clear() noexcept;

 private:
  [[nodiscard]] static Index read_index(HiCFileReader& hfs, const HiCFooter& footer);
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

#include "../../../block_reader_impl.hpp"
