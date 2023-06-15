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
#include "hictk/hic/block_cache.hpp"
#include "hictk/hic/hic_file_stream.hpp"
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
  std::shared_ptr<HiCFileStream> _hfs{};
  std::shared_ptr<BlockLRUCache> _blk_cache{};  // This should be passed in by file. Key should be
                                                // changed from size_t to {chrom1, chrom2, size_t}
  // We need the entire bin table in order to map pixels to abs bin ids
  std::shared_ptr<const BinTable> _bins{};
  Index _index{};

  BinaryBuffer _bbuffer{};
  std::vector<SerializedPixel> _tmp_buffer{};

 public:
  HiCBlockReader() = default;
  HiCBlockReader(std::shared_ptr<HiCFileStream> hfs, const Index& master_index,
                 std::shared_ptr<const BinTable> bins_,
                 std::shared_ptr<BlockLRUCache> block_cache_);

  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] const Chromosome& chrom1() const noexcept;
  [[nodiscard]] const Chromosome& chrom2() const noexcept;
  [[nodiscard]] const BinTable& bins() const noexcept;
  [[nodiscard]] const Index& index() const noexcept;

  [[nodiscard]] double sum() const noexcept;
  [[nodiscard]] double avg() const noexcept;

  [[nodiscard]] std::shared_ptr<const InteractionBlock> read(const BlockIndex& idx);

 private:
  [[nodiscard]] static Index read_index(HiCFileStream& hfs, const HiCFooter& footer);
  static void read_dispatcher_type1_block(bool i16Bin1, bool i16Bin2, bool i16Counts,
                                          std::int32_t bin1Offset, std::int32_t bin2Offset,
                                          BinaryBuffer& src,
                                          std::vector<SerializedPixel>& dest) noexcept;
  template <typename Bin1Type, typename Bin2Type, typename CountType>
  static void read_type1_block(std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer& src,
                               std::vector<SerializedPixel>& dest) noexcept;

  template <typename CountType>
  static void read_type2_block(std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer& src,
                               std::vector<SerializedPixel>& dest) noexcept;
};

}  // namespace hictk::hic::internal

#include "../../../block_reader_impl.hpp"
