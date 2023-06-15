// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cstring>
#include <memory>
#include <type_traits>
#include <utility>

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
                                      std::shared_ptr<BlockLRUCache> block_cache_)
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

inline double HiCBlockReader::avg() const noexcept {
  const auto num_bins1 = bins().subset(chrom1()).size();
  const auto num_bins2 = bins().subset(chrom2()).size();

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

inline std::shared_ptr<const InteractionBlock> HiCBlockReader::read(const BlockIndex &idx) {
  if (!idx) {
    return {nullptr};
  }

  assert(_blk_cache);
  assert(_bins);
  if (auto it = _blk_cache->find(idx.id()); it != _blk_cache->end()) {
    return it->second;
  }

  _hfs->readAndInflate(idx, _bbuffer.reset());

  const auto nRecords = static_cast<std::size_t>(_bbuffer.read<std::int32_t>());
  _tmp_buffer.resize(nRecords);

  // if (_fs->version() == 6) {
  //     readBlockOfInteractionsV6(_bbuffer, buffer);
  //     auto it =
  //         _blockCache.emplace(idx.position),
  //         InteractionBlock(buffer));
  //     return it.first->second;
  // }

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

  auto it = _blk_cache->emplace(idx.id(), InteractionBlock{idx.id(), _tmp_buffer});
  return it.first->second;
}

inline void HiCBlockReader::read_dispatcher_type1_block(
    bool i16Bin1, bool i16Bin2, bool i16Counts, std::int32_t bin1Offset, std::int32_t bin2Offset,
    BinaryBuffer &src, std::vector<SerializedPixel> &dest) noexcept {
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
      dest.push_back(SerializedPixel{bin1, bin2, counts});
    }
  }

  std::ignore = expectedNumRecords;
  assert(expectedNumRecords == dest.size());
}

template <typename CountType>
inline void HiCBlockReader::read_type2_block(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                             BinaryBuffer &src,
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

}  // namespace hictk::hic::internal
