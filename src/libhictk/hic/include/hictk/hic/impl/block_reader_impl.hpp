// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T BinaryBuffer::read() {
  static_assert(sizeof(char) == 1, "");
  assert(_i < _buffer.size());
  T x{};

  std::memcpy(static_cast<void *>(&x), _buffer.data() + _i, sizeof(T));
  _i += sizeof(T);
  return x;
}

inline std::size_t BinaryBuffer::operator()() const noexcept { return _i; }

inline std::string &BinaryBuffer::reset() noexcept {
  _buffer.clear();
  _i = 0;
  return _buffer;
}

inline HiCBlockReader::HiCBlockReader(std::shared_ptr<HiCFileReader> hfs, const Index &master_index,
                                      std::shared_ptr<const BinTable> bins_,
                                      std::shared_ptr<BlockCache> block_cache_)
    : _hfs(std::move(hfs)),
      _blk_cache(std::move(block_cache_)),
      _bins(std::move(bins_)),
      _index(master_index) {}

inline HiCBlockReader::operator bool() const noexcept { return !!_hfs; }

inline const Chromosome &HiCBlockReader::chrom1() const noexcept { return _index.chrom1(); }
inline const Chromosome &HiCBlockReader::chrom2() const noexcept { return _index.chrom2(); }

inline const BinTable &HiCBlockReader::bins() const noexcept { return *_bins; }

inline const Index &HiCBlockReader::index() const noexcept { return _index; }

inline double HiCBlockReader::sum() const noexcept { return _index.matrix_sum(); }

inline double HiCBlockReader::avg() const {
  if (_index.is_intra()) {
    throw std::domain_error(
        "HiCBlockReader::avg is not implemented for intra-chromosomal matrices");
  }

  const auto bin_size = bins().bin_size();
  const auto num_bins1 = (chrom1().size() + bin_size - 1) / bin_size;
  const auto num_bins2 = (chrom2().size() + bin_size - 1) / bin_size;

  return sum() / double(num_bins1 * num_bins2);
}

inline Index HiCBlockReader::read_index(HiCFileReader &hfs, const HiCFooter &footer) {
  if (footer.fileOffset() == -1) {
    // Footer does not exist. However, query may be valid
    return {};
  }

  return hfs.read_index(footer.fileOffset(), footer.chrom1(), footer.chrom2(), footer.unit(),
                        footer.resolution());
}

inline std::shared_ptr<const InteractionBlock> HiCBlockReader::read_v6(const Chromosome &chrom1,
                                                                       const Chromosome &chrom2,
                                                                       const BlockIndex &idx,
                                                                       bool cache_block) {
  if (!idx) {
    return {nullptr};
  }

  assert(_blk_cache);
  assert(_bins);
  auto blk = _blk_cache->find(chrom1.id(), chrom2.id(), idx.id());
  if (blk) {
    return blk;
  }

  _hfs->readAndInflate(idx, _bbuffer.reset());

  const auto nRecords = static_cast<std::size_t>(_bbuffer.read<std::int32_t>());
  _tmp_buffer.resize(nRecords);
  std::generate(_tmp_buffer.begin(), _tmp_buffer.end(), [&]() -> ThinPixel<float> {
    return {static_cast<std::uint64_t>(_bbuffer.read<std::int32_t>()),
            static_cast<std::uint64_t>(_bbuffer.read<std::int32_t>()), _bbuffer.read<float>()};
  });

  if (!cache_block) {
    return std::make_shared<const InteractionBlock>(
        InteractionBlock{idx.id(), _index.block_bin_count(), std::move(_tmp_buffer)});
  }

  return _blk_cache->emplace(
      chrom1.id(), chrom2.id(), idx.id(),
      InteractionBlock{idx.id(), _index.block_bin_count(), std::move(_tmp_buffer)});
}

inline std::shared_ptr<const InteractionBlock> HiCBlockReader::read(const Chromosome &chrom1,
                                                                    const Chromosome &chrom2,
                                                                    const BlockIndex &idx,
                                                                    bool cache_block) {
  if (_hfs->version() == 6) {
    return read_v6(chrom1, chrom2, idx, cache_block);
  }

  if (!idx) {
    return {nullptr};
  }

  assert(_blk_cache);
  assert(_bins);
  auto blk = _blk_cache->find(chrom1.id(), chrom2.id(), idx.id());
  if (blk) {
    return blk;
  }

  _hfs->readAndInflate(idx, _bbuffer.reset());

  const auto nRecords = static_cast<std::size_t>(_bbuffer.read<std::int32_t>());
  _tmp_buffer.resize(nRecords);

  const auto bin1Offset = _bbuffer.read<std::int32_t>();
  const auto bin2Offset = _bbuffer.read<std::int32_t>();

  const auto i16Counts = _bbuffer.read<char>() == 0;

  auto readUseShortBinFlag = [&]() {
    if (_hfs->version() > 8) {
      return _bbuffer.read<char>() == 0;
    }
    return true;
  };

  const auto i16Bin1 = readUseShortBinFlag();
  const auto i16Bin2 = readUseShortBinFlag();

  const auto type = static_cast<std::int8_t>(_bbuffer.read<char>());
  if (type != 1 && type != 2) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("uknown interaction type \"{}\". Supported types: 1, 2"), type));
  }

  switch (type) {
    case 1:
      read_dispatcher_type1_block(i16Bin1, i16Bin2, i16Counts, bin1Offset, bin2Offset, _bbuffer,
                                  _tmp_buffer);
      break;
    case 2:
      if (i16Counts) {
        read_type2_block<std::int16_t>(bin1Offset, bin2Offset, _bbuffer, _tmp_buffer);
        break;
      }
      read_type2_block<float>(bin1Offset, bin2Offset, _bbuffer, _tmp_buffer);
      break;
    default:
      HICTK_UNREACHABLE_CODE;
  }

  if (!cache_block) {
    return std::make_shared<const InteractionBlock>(
        InteractionBlock{idx.id(), _index.block_bin_count(), std::move(_tmp_buffer)});
  }

  return _blk_cache->emplace(
      chrom1.id(), chrom2.id(), idx.id(),
      InteractionBlock{idx.id(), _index.block_bin_count(), std::move(_tmp_buffer)});
}

inline std::size_t HiCBlockReader::read_size(const Chromosome &chrom1, const Chromosome &chrom2,
                                             const BlockIndex &idx) {
  if (!idx) {
    return 0;
  }

  assert(_blk_cache);
  assert(_bins);
  auto blk = _blk_cache->find(chrom1.id(), chrom2.id(), idx.id());
  if (blk) {
    return blk->size();
  }

  _hfs->readAndInflate(idx, _bbuffer.reset());

  return static_cast<std::size_t>(_bbuffer.read<std::int32_t>());
}

inline void HiCBlockReader::evict(const InteractionBlock &blk) {
  _blk_cache->try_erase(chrom1().id(), chrom2().id(), blk.id());
}

inline void HiCBlockReader::evict(const Chromosome &chrom1, const Chromosome &chrom2,
                                  const BlockIndex &idx) {
  _blk_cache->try_erase(chrom1.id(), chrom2.id(), idx.id());
}

inline void HiCBlockReader::clear() noexcept { _blk_cache->clear(); }

inline std::size_t HiCBlockReader::cache_size() const noexcept { return _blk_cache->size(); }

inline void HiCBlockReader::read_dispatcher_type1_block(
    bool i16Bin1, bool i16Bin2, bool i16Counts, std::int32_t bin1Offset, std::int32_t bin2Offset,
    BinaryBuffer &src, std::vector<ThinPixel<float>> &dest) noexcept {
  using BS = std::int16_t;  // Short type for bins
  using CS = std::int16_t;  // Short type for count

  using BL = std::int32_t;  // Long type for bins
  using CL = float;         // Long type for count

  if (i16Bin1 && i16Bin2 && i16Counts) {
    read_type1_block<BS, BS, CS>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (!i16Bin1 && i16Bin2 && i16Counts) {
    read_type1_block<BL, BS, CS>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (i16Bin1 && !i16Bin2 && i16Counts) {
    read_type1_block<BS, BL, CS>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (i16Bin1 && i16Bin2 && !i16Counts) {
    read_type1_block<BS, BS, CL>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (!i16Bin1 && !i16Bin2 && i16Counts) {
    read_type1_block<BL, BL, CS>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (!i16Bin1 && i16Bin2 && !i16Counts) {
    read_type1_block<BL, BS, CL>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  if (i16Bin1 && !i16Bin2 && !i16Counts) {
    read_type1_block<BS, BL, CL>(bin1Offset, bin2Offset, src, dest);
    return;
  }
  assert(!i16Bin1 && !i16Bin2 && !i16Counts);
  read_type1_block<BL, BL, CL>(bin1Offset, bin2Offset, src, dest);
}

template <typename Bin1Type, typename Bin2Type, typename CountType>
inline void HiCBlockReader::read_type1_block(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                             BinaryBuffer &src,
                                             std::vector<ThinPixel<float>> &dest) noexcept {
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
  assert(src() == expectedOffsetV7 || src() == expectedOffsetV8plus);

  const auto expectedNumRecords = dest.size();
  dest.clear();
  const auto numRows = static_cast<i32>(src.read<Bin2Type>());
  for (i32 i = 0; i < numRows; ++i) {
    const auto bin2 = bin2Offset + static_cast<i32>(src.read<Bin2Type>());

    const auto numCols = static_cast<i32>(src.read<Bin1Type>());
    for (i32 j = 0; j < numCols; ++j) {
      const auto bin1 = bin1Offset + static_cast<i32>(src.read<Bin1Type>());

      const auto counts = static_cast<f32>(src.read<CountType>());
      dest.push_back(ThinPixel<float>{static_cast<std::uint64_t>(bin1),
                                      static_cast<std::uint64_t>(bin2), counts});
    }
  }

  std::ignore = expectedNumRecords;
  assert(expectedNumRecords == dest.size());
}

template <typename CountType>
inline void HiCBlockReader::read_type2_block(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                             BinaryBuffer &src,
                                             std::vector<ThinPixel<float>> &dest) noexcept {
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
    const auto bin1 = static_cast<std::uint64_t>(bin1Offset + col);
    const auto bin2 = static_cast<std::uint64_t>(bin2Offset + row);

    dest.emplace_back(ThinPixel<float>{bin1, bin2, static_cast<f32>(count)});
  }
}

}  // namespace hictk::hic::internal
