// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/block_reader.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/hic/interaction_block.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

HiCBlockReader::HiCBlockReader(std::shared_ptr<HiCFileReader> hfs, Index master_index,
                               std::shared_ptr<const BinTable> bins_,
                               std::shared_ptr<BlockCache> block_cache_) noexcept
    : _hfs(std::move(hfs)),
      _blk_cache(std::move(block_cache_)),
      _bins(std::move(bins_)),
      _index(std::move(master_index)) {}

HiCBlockReader::operator bool() const noexcept { return !!_hfs; }

const Chromosome &HiCBlockReader::chrom1() const noexcept { return _index.chrom1(); }
const Chromosome &HiCBlockReader::chrom2() const noexcept { return _index.chrom2(); }

const BinTable &HiCBlockReader::bins() const noexcept {
  assert(_bins);
  return *_bins;
}

std::shared_ptr<const BinTable> HiCBlockReader::bins_ptr() const noexcept { return _bins; }

const Index &HiCBlockReader::index() const noexcept { return _index; }

double HiCBlockReader::sum() const noexcept { return _index.matrix_sum(); }

double HiCBlockReader::avg() const {
  if (_index.is_intra()) {
    throw std::domain_error(
        "HiCBlockReader::avg is not implemented for intra-chromosomal matrices");
  }

  const auto bin_size = bins().resolution();
  // We round down for two reasons:
  // - to be consistent with straw
  // - because the last bin is usually smaller than the bin_size
  const auto num_bins1 = chrom1().size() / bin_size;
  const auto num_bins2 = chrom2().size() / bin_size;

  return sum() / static_cast<double>(num_bins1 * num_bins2);
}

std::shared_ptr<const InteractionBlock> HiCBlockReader::read_v6(const Chromosome &chrom1,
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

std::shared_ptr<const InteractionBlock> HiCBlockReader::read(const Chromosome &chrom1,
                                                             const Chromosome &chrom2,
                                                             const BlockIndex &idx,
                                                             bool cache_block) {
  if (_hfs->version() == 6) {  // NOLINT(*-avoid-magic-numbers)
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
    if (_hfs->version() > 8) {  // NOLINT(*-avoid-magic-numbers)
      return _bbuffer.read<char>() == 0;
    }
    return true;
  };

  const auto i16Bin1 = readUseShortBinFlag();
  const auto i16Bin2 = readUseShortBinFlag();

  const auto type = static_cast<std::int8_t>(_bbuffer.read<char>());
  if (type != 1 && type != 2) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unknown interaction type \"{}\". Supported types: 1, 2"), type));
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
      unreachable_code();
  }

  if (!cache_block) {
    return std::make_shared<const InteractionBlock>(
        InteractionBlock{idx.id(), _index.block_bin_count(), std::move(_tmp_buffer)});
  }

  return _blk_cache->emplace(
      chrom1.id(), chrom2.id(), idx.id(),
      InteractionBlock{idx.id(), _index.block_bin_count(), std::move(_tmp_buffer)});
}

std::size_t HiCBlockReader::read_size(const Chromosome &chrom1, const Chromosome &chrom2,
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

void HiCBlockReader::evict(const InteractionBlock &blk) {
  _blk_cache->try_erase(chrom1().id(), chrom2().id(), blk.id());
}

void HiCBlockReader::evict(const Chromosome &chrom1, const Chromosome &chrom2,
                           const BlockIndex &idx) {
  _blk_cache->try_erase(chrom1.id(), chrom2.id(), idx.id());
}

void HiCBlockReader::clear() noexcept { _blk_cache->clear(); }

std::size_t HiCBlockReader::cache_size() const noexcept { return _blk_cache->size(); }

void HiCBlockReader::read_dispatcher_type1_block(bool i16Bin1, bool i16Bin2, bool i16Counts,
                                                 std::int32_t bin1Offset, std::int32_t bin2Offset,
                                                 BinaryBuffer &src,
                                                 std::vector<ThinPixel<float>> &dest) noexcept {
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

}  // namespace hictk::hic::internal
