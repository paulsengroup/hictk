// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <type_traits>
#include <vector>

#include "hictk/hic/cache.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/hic_file_stream.hpp"

namespace hictk::internal {

class MatrixSelector {
  struct BinaryBuffer {
    std::string buffer{};
    std::size_t i{};

    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    T read();
  };

  std::shared_ptr<HiCFileStream> _fs;
  std::shared_ptr<const HiCFooter> _footer;
  BlockMap _blockMap{};
  BlockLRUCache _blockCache{};
  std::set<std::size_t> _blockNumberBuff{};
  std::vector<contactRecord> _contactRecordBuff{};
  BinaryBuffer _buffer{};

 public:
  MatrixSelector() = delete;
  MatrixSelector(std::shared_ptr<HiCFileStream> fs, std::shared_ptr<const HiCFooter> footer,
                 std::size_t block_cache_capacity);

  [[nodiscard]] const chromosome &chrom1() const noexcept;
  [[nodiscard]] const chromosome &chrom2() const noexcept;

  [[nodiscard]] std::int64_t resolution() const noexcept;
  [[nodiscard]] MatrixType matrix_type() const noexcept;
  [[nodiscard]] NormalizationMethod normalizationMethod() const noexcept;
  [[nodiscard]] MatrixUnit matrixUnit() const noexcept;

  [[nodiscard]] std::int64_t numBins1() const noexcept;
  [[nodiscard]] std::int64_t numBins2() const noexcept;

  [[nodiscard]] bool isIntra() const noexcept;
  [[nodiscard]] bool isInter() const noexcept;

  [[nodiscard]] const std::vector<double> &chrom1Norm() const noexcept;
  [[nodiscard]] const std::vector<double> &chrom2Norm() const noexcept;

  [[nodiscard]] inline double avgCount() const;

  void fetch(std::vector<contactRecord> &buffer, bool sorted = false);
  void fetch(const std::string &coord, std::vector<contactRecord> &buffer, bool sorted = false);
  void fetch(std::int64_t start, std::int64_t end, std::vector<contactRecord> &buffer,
             bool sorted = false);

  void fetch(const std::string &coord1, const std::string &coord2,
             std::vector<contactRecord> &buffer, bool sorted = false);
  void fetch(std::int64_t start1, std::int64_t end1, std::int64_t start2, std::int64_t end2,
             std::vector<contactRecord> &buffer, bool sorted = false);

  void clearBlockCache() noexcept;
  [[nodiscard]] constexpr double blockCacheHitRate() const noexcept;
  [[nodiscard]] std::size_t blockCacheSize() const noexcept;
  [[nodiscard]] constexpr std::size_t blockCacheSizeBytes() const noexcept;
  [[nodiscard]] constexpr std::size_t blockCacheHits() const noexcept;
  [[nodiscard]] constexpr std::size_t blockCacheMisses() const noexcept;

 private:
  [[nodiscard]] static BlockMap readBlockMap(HiCFileStream &fs, const HiCFooter &footer);

  void readBlockNumbers(std::int64_t bin1, std::int64_t bin2, std::int64_t bin3, std::int64_t bin4,
                        std::set<std::size_t> &buffer) const;
  void readBlockNumbersV9Intra(std::int64_t bin1, std::int64_t bin2, std::int64_t bin3,
                               std::int64_t bin4, std::set<std::size_t> &buffer) const;
  [[nodiscard]] std::shared_ptr<InteractionBlock> readBlockOfInteractions(
      indexEntry idx, std::vector<contactRecord> &buffer);
  [[nodiscard]] contactRecord processInteraction(contactRecord record);
  // static void readBlockOfInteractionsV6(BinaryBuffer &src, std::vector<contactRecord> &dest);

  static void readBlockOfInteractionsType1Dispatcher(bool i16Bin1, bool i16Bin2, bool i16Counts,
                                                     std::int32_t bin1Offset,
                                                     std::int32_t bin2Offset, BinaryBuffer &src,
                                                     std::vector<contactRecord> &dest) noexcept;
  template <typename Bin1Type, typename Bin2Type, typename CountType>
  static void readBlockOfInteractionsType1(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                           BinaryBuffer &src,
                                           std::vector<contactRecord> &dest) noexcept;

  template <typename CountType>
  static void readBlockOfInteractionsType2(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                           BinaryBuffer &src,
                                           std::vector<contactRecord> &dest) noexcept;
};

}  // namespace hictk::internal

#include "../../../hic_matrix_selector_impl.hpp"
