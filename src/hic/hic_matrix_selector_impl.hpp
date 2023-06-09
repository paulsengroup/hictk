// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ios>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "hictk/hic/common.hpp"

namespace hictk::internal {

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T MatrixSelector::BinaryBuffer::read() {
  static_assert(sizeof(char) == 1, "");
  assert(i < buffer.size());
  T x{};

  std::memcpy(static_cast<void *>(&x), buffer.data() + i, sizeof(T));
  i += sizeof(T);
  return x;
}

inline MatrixSelector::MatrixSelector(std::shared_ptr<HiCFileStream> fs,
                                      std::shared_ptr<const HiCFooter> footer,
                                      std::size_t block_cache_capacity)
    : _fs(std::move(fs)),
      _footer(std::move(footer)),
      _bins(_fs->header().chromosomes, _footer->resolution()),
      _blockMap(readBlockMap(*_fs, *_footer)),
      _blockCache(block_cache_capacity) {
  assert(_footer);
}

inline const Chromosome &MatrixSelector::chrom1() const noexcept { return _footer->chrom1(); }

inline const Chromosome &MatrixSelector::chrom2() const noexcept { return _footer->chrom2(); }

inline std::uint32_t MatrixSelector::resolution() const noexcept { return _footer->resolution(); }

inline MatrixType MatrixSelector::matrix_type() const noexcept { return _footer->matrix_type(); }

inline NormalizationMethod MatrixSelector::normalizationMethod() const noexcept {
  return _footer->normalization();
}

inline MatrixUnit MatrixSelector::matrixUnit() const noexcept { return _footer->unit(); }

inline std::int64_t MatrixSelector::numBins1() const noexcept {
  return (chrom1().size() + resolution() - 1) / resolution();
}

inline std::int64_t MatrixSelector::numBins2() const noexcept {
  return (chrom2().size() + resolution() - 1) / resolution();
}

inline bool MatrixSelector::isIntra() const noexcept { return chrom1() == chrom2(); }

inline bool MatrixSelector::isInter() const noexcept { return !isIntra(); }

inline const std::vector<double> &MatrixSelector::chrom1Norm() const noexcept {
  return _footer->c1Norm();
}

inline const std::vector<double> &MatrixSelector::chrom2Norm() const noexcept {
  return _footer->c2Norm();
}

inline double MatrixSelector::avgCount() const {
  if (isInter()) {
    return _blockMap.sumCount / static_cast<double>(numBins1() * numBins2());
  }
  throw std::domain_error(
      "MatrixSelector::avgCount is not implemented for intra-chromosomal matrices");
}

inline void MatrixSelector::fetch(std::vector<Pixel<float>> &buffer, bool sorted) {
  return fetch(0, chrom1().size(), 0, chrom2().size(), buffer, sorted);
}

inline void MatrixSelector::fetch(std::int64_t start, std::int64_t end,
                                  std::vector<Pixel<float>> &buffer, bool sorted) {
  return fetch(start, end, start, end, buffer, sorted);
}

inline void MatrixSelector::fetch(std::int64_t start1, std::int64_t end1, std::int64_t start2,
                                  std::int64_t end2, std::vector<Pixel<float>> &buffer,
                                  bool sorted) {
  buffer.clear();
  if (start1 > end1) {
    throw std::invalid_argument(fmt::format(FMT_STRING("start1 > end1: {} > {}"), start1, end1));
  }
  if (start2 > end2) {
    throw std::invalid_argument(fmt::format(FMT_STRING("start2 > end2: {} > {}"), start2, end2));
  }

  if (start1 < 0 || end1 > chrom1().size()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("query extends past chromosome {}: interval {}-{} lies outside of 0-{}"),
        chrom1().name(), start1, end1, chrom1().size()));
  }

  if (start2 < 0 || end2 > chrom2().size()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("query extends past chromosome {}: interval {}-{} lies outside of 0-{}"),
        chrom2().name(), start2, end2, chrom2().size()));
  }

  // Query is valid but returns no pixels
  if (this->_footer->fileOffset() == -1) {
    assert(this->_blockMap.blocks.empty());
    return;
  }

  const auto is_intra = isIntra();
  if (is_intra && start1 > start2) {
    std::swap(start1, start2);
    std::swap(end1, end2);
  }

  const auto bin1 = start1 / resolution();
  const auto bin2 = (end1 + resolution() - 1) / resolution();
  const auto bin3 = start2 / resolution();
  const auto bin4 = (end2 + resolution() - 1) / resolution();

  if (_fs->version() > 8 && isIntra()) {
    readBlockNumbersV9Intra(bin1, bin2, bin3, bin4, _blockNumberBuff);
  } else {
    readBlockNumbers(bin1, bin2, bin3, bin4, _blockNumberBuff);
  }

  std::size_t empty_blocks = 0;
  for (auto blockNumber : _blockNumberBuff) {
    const auto block = readBlockOfInteractions(_blockMap.blocks[blockNumber], _contactRecordBuff);
    if (!block) {
      empty_blocks++;
      continue;
    }

    // Obs we use open-closed interval instead of open-open like is done in straw
    for (const auto &[b1, row] : block->find_overlap(bin1, bin2)) {
      if (b1 >= bin2) {
        // We're past the last row overlapping the query
        break;
      }
      for (const auto &tp : row) {
        const auto &b2 = tp.bin2_id;
        if (b1 < bin1 || b2 < bin3) {
          // We're upstream of the first column overlapping the query (if any)
          continue;
        }

        if (b2 >= bin4) {
          // We're past the last column overlapping the query for the current row
          break;
        }

        auto record = processInteraction(SerializedPixel{b1, b2, tp.count});
        if (std::isfinite(record.count)) {
          buffer.emplace_back(
              PixelCoordinates{
                  _bins.at(_footer->chrom1(), static_cast<std::uint32_t>(record.bin1_id)),
                  _bins.at(_footer->chrom2(), static_cast<std::uint32_t>(record.bin2_id))},
              record.count);
        }
      }
    };
  }
  if (sorted && _blockNumberBuff.size() - empty_blocks > 1) {
    // Only interactions from the same block are guaranteed to already be sorted
    std::sort(buffer.begin(), buffer.end());
  }
}

inline SerializedPixel MatrixSelector::processInteraction(SerializedPixel record) {
  const auto &c1Norm = _footer->c1Norm();
  const auto &c2Norm = _footer->c2Norm();
  const auto &expected = _footer->expectedValues();

  assert(isInter() || record.bin1_id <= record.bin2_id);

  const auto skipNormalization =
      normalizationMethod() == NormalizationMethod::NONE || matrix_type() == MatrixType::expected;

  if (!skipNormalization) {
    const auto bin1 = static_cast<std::size_t>(record.bin1_id);
    const auto bin2 = static_cast<std::size_t>(record.bin2_id);
    assert(bin1 < c1Norm.size());
    assert(bin2 < c2Norm.size());
    record.count /= static_cast<float>(c1Norm[bin1] * c2Norm[bin2]);
  }

  record.bin1_id *= resolution();
  record.bin2_id *= resolution();

  if (matrix_type() == MatrixType::observed) {
    return record;
  }

  const auto expectedCount = [&]() {
    if (isInter()) {
      return float(avgCount());
    }

    const auto i = static_cast<std::size_t>((record.bin2_id - record.bin1_id) / resolution());
    assert(i < expected.size());
    return float(expected[i]);
  }();

  if (matrix_type() == MatrixType::expected) {
    record.count = expectedCount;
    return record;
  }

  assert(matrix_type() == MatrixType::oe);
  record.count /= expectedCount;

  return record;
}

/*
inline void MatrixSelector::readBlockOfInteractionsV6(BinaryBuffer &src,
                                                      std::vector<SerializedPixel> &dest) {
    assert(src.i == sizeof(std::int32_t));

    constexpr auto recordSize = sizeof(std::int32_t) + sizeof(std::int32_t) + sizeof(float);

    const auto srcSize = sizeof(std::int32_t) + (src.buffer.size() * sizeof(char));
    const auto destSize = recordSize * dest.size();

    if (srcSize != destSize) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("binary buffer appears to be corrupted: expected {}B, found {}B"), destSize,
            srcSize));
    }

    std::generate(dest.begin(), dest.end(), [&]() {
        // clang-format off
        return SerializedPixel{src.read<std::int32_t>(),
                             src.read<std::int32_t>(),
                             src.read<float>()};
        // clang-format on
    });
    return;
}
*/

inline std::shared_ptr<InteractionBlock> MatrixSelector::readBlockOfInteractions(
    indexEntry idx, std::vector<SerializedPixel> &buffer) {
  buffer.clear();
  if (idx.size <= 0) {
    return {nullptr};
  }

  if (auto it = _blockCache.find(static_cast<std::size_t>(idx.position)); it != _blockCache.end()) {
    return it->second;
  }

  _fs->readAndInflate(idx, _buffer.buffer);
  _buffer.i = 0;

  const auto nRecords = static_cast<std::size_t>(_buffer.read<std::int32_t>());
  buffer.resize(nRecords);

  // if (_fs->version() == 6) {
  //     readBlockOfInteractionsV6(_buffer, buffer);
  //     auto it =
  //         _blockCache.emplace(static_cast<std::size_t>(idx.position),
  //         InteractionBlock(buffer));
  //     return it.first->second;
  // }

  const auto bin1Offset = _buffer.read<std::int32_t>();
  const auto bin2Offset = _buffer.read<std::int32_t>();

  const auto i16Counts = _buffer.read<char>() == 0;

  auto readUseShortBinFlag = [&]() {
    if (_fs->version() > 8) {
      return _buffer.read<char>() == 0;
    }
    return true;
  };

  const auto i16Bin1 = readUseShortBinFlag();
  const auto i16Bin2 = readUseShortBinFlag();

  const auto type = static_cast<std::int8_t>(_buffer.read<char>());
  if (type != 1 && type != 2) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("uknown interaction type \"{}\". Supported types: 1, 2"), type));
  }

  switch (type) {
    case 1:
      readBlockOfInteractionsType1Dispatcher(i16Bin1, i16Bin2, i16Counts, bin1Offset, bin2Offset,
                                             _buffer, buffer);
      break;
    case 2:
      if (i16Counts) {
        readBlockOfInteractionsType2<std::int16_t>(bin1Offset, bin2Offset, _buffer, buffer);
        break;
      }
      readBlockOfInteractionsType2<float>(bin1Offset, bin2Offset, _buffer, buffer);
      break;
    default:
      assert(false);
      std::abort();
  }

  auto it = _blockCache.emplace(static_cast<std::size_t>(idx.position), InteractionBlock(buffer));
  return it.first->second;
}

inline void MatrixSelector::readBlockOfInteractionsType1Dispatcher(
    bool i16Bin1, bool i16Bin2, bool i16Counts, std::int32_t bin1Offset, std::int32_t bin2Offset,
    BinaryBuffer &src, std::vector<SerializedPixel> &dest) noexcept {
  using BS = std::int16_t;  // Short type for bins
  using CS = std::int16_t;  // Short type for count

  using BL = std::int32_t;  // Long type for bins
  using CL = float;         // Long type for count

  if (i16Bin1 && i16Bin2 && i16Counts) {
    readBlockOfInteractionsType1<BS, BS, CS>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (!i16Bin1 && i16Bin2 && i16Counts) {
    readBlockOfInteractionsType1<BL, BS, CS>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (i16Bin1 && !i16Bin2 && i16Counts) {
    readBlockOfInteractionsType1<BS, BL, CS>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (i16Bin1 && i16Bin2 && !i16Counts) {
    readBlockOfInteractionsType1<BS, BS, CL>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (!i16Bin1 && !i16Bin2 && i16Counts) {
    readBlockOfInteractionsType1<BL, BL, CS>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (!i16Bin1 && i16Bin2 && !i16Counts) {
    readBlockOfInteractionsType1<BL, BS, CL>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (i16Bin1 && !i16Bin2 && !i16Counts) {
    readBlockOfInteractionsType1<BS, BL, CL>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  assert(!i16Bin1 && !i16Bin2 && !i16Counts);
  readBlockOfInteractionsType1<BL, BL, CL>(bin1Offset, bin2Offset, src, dest);
}

template <typename Bin1Type, typename Bin2Type, typename CountType>
inline void MatrixSelector::readBlockOfInteractionsType1(
    std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer &src,
    std::vector<SerializedPixel> &dest) noexcept {
  using i16 = std::int16_t;
  using i32 = std::int32_t;
  using f32 = float;
  static_assert(std::is_same<i16, Bin1Type>::value || std::is_same<i32, Bin1Type>::value, "");
  static_assert(std::is_same<i16, Bin2Type>::value || std::is_same<i32, Bin2Type>::value, "");
  static_assert(std::is_same<i16, CountType>::value || std::is_same<f32, CountType>::value, "");

  constexpr auto expectedOffsetV7 = (3 * sizeof(i32)) + (2 * sizeof(char));
  constexpr auto expectedOffsetV8plus = expectedOffsetV7 + (2 * sizeof(char));
  std::ignore = expectedOffsetV7;
  std::ignore = expectedOffsetV8plus;
  assert(src.i == expectedOffsetV7 || src.i == expectedOffsetV8plus);

  const auto expectedNumRecords = dest.size();
  dest.clear();
  const auto numRows = static_cast<i32>(src.read<Bin2Type>());
  for (i32 i = 0; i < numRows; ++i) {
    const auto bin2 = bin2Offset + static_cast<i32>(src.read<Bin2Type>());

    const auto numCols = static_cast<i32>(src.read<Bin1Type>());
    for (i32 j = 0; j < numCols; ++j) {
      const auto bin1 = bin1Offset + static_cast<i32>(src.read<Bin1Type>());

      const auto counts = static_cast<f32>(src.read<CountType>());
      dest.push_back(SerializedPixel{bin1, bin2, counts});
    }
  }

  std::ignore = expectedNumRecords;
  assert(expectedNumRecords == dest.size());
}

template <typename CountType>
inline void MatrixSelector::readBlockOfInteractionsType2(
    std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer &src,
    std::vector<SerializedPixel> &dest) noexcept {
  using i16 = std::int16_t;
  using i32 = std::int32_t;
  using f32 = float;
  static_assert(std::is_same<i16, CountType>::value || std::is_same<f32, CountType>::value, "");

  const auto nPts = src.read<i32>();
  const auto w = static_cast<i32>(src.read<i16>());

  constexpr auto i16Sentinel = (std::numeric_limits<i16>::lowest)();
  constexpr auto i16Counts = std::is_same<i16, CountType>::value;

  auto isValid = [&](CountType n) {
    return (i16Counts && static_cast<i16>(n) != i16Sentinel) ||
           (!i16Counts && !std::isnan(static_cast<f32>(n)));
  };

  dest.reserve(static_cast<std::size_t>(nPts));
  dest.clear();
  for (i32 i = 0; i < nPts; ++i) {
    const auto count = src.read<CountType>();
    if (!isValid(count)) {
      continue;
    }
    const auto row = i / w;
    const auto col = i - row * w;
    const auto bin1 = bin1Offset + col;
    const auto bin2 = bin2Offset + row;

    dest.emplace_back(SerializedPixel{bin1, bin2, static_cast<f32>(count)});
  }
}

inline BlockMap MatrixSelector::readBlockMap(HiCFileStream &fs, const HiCFooter &footer) {
  if (footer.fileOffset() == -1) {
    // Footer does not exist. However, query validity is assessed elswehere
    return {};
  }

  BlockMap buffer{};
  fs.readBlockMap(footer.fileOffset(), footer.chrom1(), footer.chrom2(), footer.unit(),
                  footer.resolution(), buffer);
  return buffer;
}

inline void MatrixSelector::readBlockNumbers(std::int64_t bin1, std::int64_t bin2,
                                             std::int64_t bin3, std::int64_t bin4,
                                             std::set<std::size_t> &buffer) const {
  const auto blockBinCount = _blockMap.blockBinCount;
  const auto blockColumnCount = _blockMap.blockColumnCount;

  const auto col1 = bin1 / blockBinCount;
  const auto col2 = (bin2 + 1) / blockBinCount;
  const auto row1 = bin3 / blockBinCount;
  const auto row2 = (bin4 + 1) / blockBinCount;

  // check region part that overlaps with lower left triangle but only if intrachromosomal
  const auto checkLowerLeftTri = isIntra();
  buffer.clear();
  // first check the upper triangular matrix_type
  for (auto row = row1; row <= row2; ++row) {
    for (auto col = col1; col <= col2; ++col) {
      buffer.insert(static_cast<std::size_t>(row * blockColumnCount + col));
      if (checkLowerLeftTri) {
        buffer.insert(static_cast<std::size_t>(col * blockColumnCount + row));
      }
    }
  }
}

inline void MatrixSelector::readBlockNumbersV9Intra(std::int64_t bin1, std::int64_t bin2,
                                                    std::int64_t bin3, std::int64_t bin4,
                                                    std::set<std::size_t> &buffer) const {
  const auto blockBinCount = _blockMap.blockBinCount;
  const auto blockColumnCount = _blockMap.blockColumnCount;

  const auto translatedLowerPAD = (bin1 + bin3) / 2 / blockBinCount;
  const auto translatedHigherPAD = (bin2 + bin4) / 2 / blockBinCount + 1;
  const auto translatedNearerDepth = static_cast<std::int64_t>(
      std::log2(1.0 + double(std::abs(bin1 - bin4)) / std::sqrt(2.0) / blockBinCount));
  const auto translatedFurtherDepth = static_cast<std::int64_t>(
      std::log2(1.0 + double(std::abs(bin2 - bin3)) / std::sqrt(2.0) / blockBinCount));

  // because code above assumes above diagonal; but we could be below diagonal
  const auto nearerDepth = [&]() -> std::int64_t {
    if ((bin1 > bin4 && bin2 < bin3) || (bin2 > bin3 && bin1 < bin4)) {
      return 0;
    }
    return std::min(translatedNearerDepth, translatedFurtherDepth);
  }();

  // +1; integer divide rounds down
  const auto furtherDepth = std::max(translatedNearerDepth, translatedFurtherDepth) + 1;

  buffer.clear();
  for (auto depth = nearerDepth; depth <= furtherDepth; ++depth) {
    for (auto pad = translatedLowerPAD; pad <= translatedHigherPAD; ++pad) {
      buffer.insert(static_cast<std::size_t>(depth * blockColumnCount + pad));
    }
  }
}

inline void MatrixSelector::clearBlockCache() noexcept { _blockCache.reset(); }
constexpr double MatrixSelector::blockCacheHitRate() const noexcept {
  return _blockCache.hit_rate();
}
inline std::size_t MatrixSelector::blockCacheSize() const noexcept { return _blockCache.size(); }
constexpr std::size_t MatrixSelector::blockCacheSizeBytes() const noexcept {
  return _blockCache.size_in_bytes();
}
constexpr std::size_t MatrixSelector::blockCacheHits() const noexcept { return _blockCache.hits(); }
constexpr std::size_t MatrixSelector::blockCacheMisses() const noexcept {
  return _blockCache.misses();
}

}  // namespace hictk::internal
