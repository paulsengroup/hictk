// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>
#include <zstd.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/hic/binary_buffer.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

template <typename N>
inline void MatrixInteractionBlockFlat<N>::emplace_back(Pixel<N> &&p) {
  emplace_back(p.to_thin());
}

template <typename N>
inline void MatrixInteractionBlockFlat<N>::emplace_back(ThinPixel<N> &&p) {
  bin1_ids.push_back(p.bin1_id);
  bin2_ids.push_back(p.bin2_id);
  counts.push_back(p.count);
}

template <typename N>
inline std::size_t MatrixInteractionBlockFlat<N>::size() const noexcept {
  return bin1_ids.size();
}

template <typename N>
inline std::string MatrixInteractionBlockFlat<N>::serialize(BinaryBuffer &buffer,
                                                            ZSTD_CCtx_s &compressor,
                                                            std::string &compression_buffer,
                                                            int compression_lvl, bool clear) const {
  if (size() == 0) {
    return "";
  }

  if (clear) {
    buffer.clear();
  }

  buffer.write(bin1_ids);
  buffer.write(bin2_ids);
  buffer.write(counts);

  const auto buff_size = ZSTD_compressBound(buffer.get().size() * sizeof(char));
  compression_buffer.resize(buff_size);

  std::size_t compressed_size = ZSTD_compressCCtx(
      &compressor, reinterpret_cast<void *>(compression_buffer.data()),
      compression_buffer.size() * sizeof(char), reinterpret_cast<const void *>(buffer.get().data()),
      buffer.get().size() * sizeof(char), compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  compression_buffer.resize(compressed_size);

  buffer.clear();
  buffer.write(size());
  buffer.write(compression_buffer, false);

  return buffer.get();
}

template <typename N>
[[nodiscard]] std::vector<ThinPixel<N>> MatrixInteractionBlockFlat<N>::deserialize(
    BinaryBuffer &buffer, ZSTD_DCtx_s &decompressor, std::string &decompression_buffer) {
  const auto size_ = buffer.read<std::size_t>();
  std::vector<ThinPixel<N>> pixels(size_);

  const auto decompressed_size =
      size_ * (sizeof(std::uint64_t) + sizeof(std::uint64_t) + sizeof(N));
  decompression_buffer.resize(decompressed_size);

  const auto compressed_buffer = std::string_view{buffer.get()}.substr(sizeof(size_));

  const auto status = ZSTD_decompressDCtx(
      &decompressor, decompression_buffer.data(), decompression_buffer.size() * sizeof(char),
      compressed_buffer.data(), compressed_buffer.size() * sizeof(char));

  if (ZSTD_isError(status)) {
    throw std::runtime_error(ZSTD_getErrorName(status));
  }
  buffer.clear();
  buffer.write(decompression_buffer);

  for (auto &p : pixels) {
    p.bin1_id = buffer.read<std::uint64_t>();
  }
  for (auto &p : pixels) {
    p.bin2_id = buffer.read<std::uint64_t>();
  }
  for (auto &p : pixels) {
    p.count = buffer.read<N>();
  }

  return pixels;
}

inline bool HiCInteractionToBlockMapper::BlockID::operator<(const BlockID &other) const noexcept {
  if (chrom1_id != other.chrom1_id) {
    return chrom1_id < other.chrom1_id;
  }
  if (chrom2_id != other.chrom2_id) {
    return chrom2_id < other.chrom2_id;
  }
  return bid < other.bid;
}

inline HiCInteractionToBlockMapper::HiCInteractionToBlockMapper(
    std::filesystem::path path, std::shared_ptr<const BinTable> bins, int compression_lvl)
    : _path(std::move(path)),
      _fs(filestream::FileStream::create(_path)),
      _bin_table(std::move(bins)),
      _compression_lvl(compression_lvl),
      _zstd_cctx(ZSTD_createCCtx()),
      _zstd_dctx(ZSTD_createDCtx()) {
  init_block_mappers();
}

inline const Reference &HiCInteractionToBlockMapper::chromosomes() const noexcept {
  return _bin_table->chromosomes();
}

template <typename PixelIt, typename>
inline void HiCInteractionToBlockMapper::append_pixels(PixelIt first_pixel, PixelIt last_pixel,
                                                       std::size_t chunk_size) {
  while (first_pixel != last_pixel) {
    if (++_pixels_processed >= chunk_size) {
      write_blocks();
      _blocks.clear();
      _pixels_processed = 0;
    }

    auto p = *first_pixel;
    auto bid = map(p);

    const auto &chrom1 = chromosomes().at(bid.chrom1_id);
    const auto &chrom2 = chromosomes().at(bid.chrom2_id);
    const auto chrom_pair = std::make_pair(chrom1, chrom2);

    auto match = _blocks.find(bid);
    if (match != _blocks.end()) {
      _pixel_sums.at(chrom_pair) += p.count;
      match->second.emplace_back(std::move(p));
    } else {
      _pixel_sums.emplace(chrom_pair, p.count);
      auto [it, _] = _blocks.emplace(std::move(bid), MatrixInteractionBlockFlat<float>{});
      it->second.emplace_back(std::move(p));
    }
    ++first_pixel;
  }
}

inline auto HiCInteractionToBlockMapper::block_index() const noexcept -> const BlockIndexMap& {
  return _block_index;
}

inline auto HiCInteractionToBlockMapper::chromosome_index() const noexcept -> const ChromosomeIndexMap& {
  return _chromosome_index;
}

inline auto HiCInteractionToBlockMapper::merge_blocks(
    const HiCInteractionToBlockMapper::BlockID &bid) -> MatrixInteractionBlock<float> {
  MatrixInteractionBlock<float> blk{};
  for (const auto &[pos, size] : _block_index.at(bid)) {
    _fs.seekg(static_cast<std::streamoff>(pos));
    _fs.read(_bbuffer.reset(), size);
    for (const auto &pixel : MatrixInteractionBlockFlat<float>::deserialize(_bbuffer, *_zstd_dctx,
                                                                            _compression_buffer)) {
      blk.emplace_back(Pixel<float>(*_bin_table, pixel));
    }
  }
  blk.finalize();
  return blk;
}

inline float HiCInteractionToBlockMapper::pixel_sum(const Chromosome &chrom1,
                                                    const Chromosome &chrom2) {
  return _pixel_sums.at(std::make_pair(chrom1, chrom2));
}

inline void HiCInteractionToBlockMapper::finalize() {
  write_blocks();
  index_chromosomes();
  _blocks.clear();
  _pixels_processed = 0;
  _fs.flush();
}

inline void HiCInteractionToBlockMapper::init_block_mappers() {
  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);

      const auto num_bins = compute_num_bins(chrom1, chrom2, _bin_table->bin_size());
      const auto num_columns = compute_block_column_count(
          chrom1, chrom2, _bin_table->bin_size(),
          chrom1 == chrom2 ? DEFAULT_INTRA_CUTOFF : DEFAULT_INTER_CUTOFF);
      const auto num_rows = num_bins / num_columns + 1;

      if (chrom1 == chrom2) {
        _mappers_intra.emplace(chrom1, BlockMapperIntra{num_rows, num_columns});
      } else {
        _mappers_inter.emplace(std::make_pair(chrom1, chrom2),
                               BlockMapperInter{num_rows, num_columns});
      }
    }
  }
}

inline std::pair<std::uint64_t, std::uint32_t> HiCInteractionToBlockMapper::write_block(
    const MatrixInteractionBlockFlat<float> &blk) {
  const auto offset = _fs.tellp();
  _fs.write(blk.serialize(_bbuffer, *_zstd_cctx, _compression_buffer, _compression_lvl));
  const auto size = _fs.tellp() - offset;
  return std::make_pair(offset, static_cast<std::uint32_t>(size));
}

template <typename N>
inline auto HiCInteractionToBlockMapper::map(const ThinPixel<N> &p) const -> BlockID {
  return map(Pixel<N>(*_bin_table, p));
}

template <typename N>
inline auto HiCInteractionToBlockMapper::map(const Pixel<N> &p) const -> BlockID {
  const auto &bin1 = p.coords.bin1;
  const auto &bin2 = p.coords.bin2;

  const auto &chrom1 = bin1.chrom();
  const auto &chrom2 = bin2.chrom();

  const auto bin1_id = bin1.rel_id();
  const auto bin2_id = bin2.rel_id();

  const auto block_id = p.coords.is_intra()
                            ? _mappers_intra.at(chrom1)(bin1_id, bin2_id)
                            : _mappers_inter.at(std::make_pair(chrom1, chrom2))(bin1_id, bin2_id);

  return {chrom1.id(), chrom2.id(), block_id};
}

inline void HiCInteractionToBlockMapper::write_blocks() {
  for (auto &[bid, blk] : _blocks) {
    const auto [offset, size] = write_block(blk);

    auto match = _block_index.find(bid);
    if (match != _block_index.end()) {
      match->second.emplace_back(BlockIndex{offset, size});
    } else {
      _block_index.emplace(bid, std::vector<BlockIndex>{{offset, size}});
    }
  }
}

inline void HiCInteractionToBlockMapper::index_chromosomes() {
  for (const auto &[bid, _] : _block_index) {
    const auto key = std::make_pair(_bin_table->chromosomes().at(bid.chrom1_id),
                                    _bin_table->chromosomes().at(bid.chrom2_id));
    auto match = _chromosome_index.find(key);
    if (match != _chromosome_index.end()) {
      match->second.emplace_back(bid);
    } else {
      _chromosome_index.emplace(key, std::vector<BlockID>{bid});
    }
  }

  for (auto &[_, v] : _chromosome_index) {
    std::sort(v.begin(), v.end());
  }
}

inline std::size_t HiCInteractionToBlockMapper::compute_block_column_count(
    const Chromosome &chrom1, const Chromosome &chrom2, std::uint32_t bin_size,
    std::uint32_t cutoff, std::size_t block_capacity) {
  const auto num_bins = compute_num_bins(chrom1, chrom2, bin_size);
  auto num_columns = num_bins / block_capacity + 1;
  if (bin_size < cutoff) {
    const auto genome_size = num_bins * bin_size;
    num_columns = genome_size / (block_capacity * cutoff);
  }

  const auto max_sqrt =
      static_cast<std::size_t>(std::sqrt(std::numeric_limits<std::int32_t>::max()));
  return std::clamp(num_columns, std::size_t(1), max_sqrt - 1);
}

inline std::size_t HiCInteractionToBlockMapper::compute_num_bins(const Chromosome &chrom1,
                                                                 const Chromosome &chrom2,
                                                                 std::size_t bin_size) {
  const auto max_size = std::max(chrom1.size(), chrom2.size());
  return (max_size + bin_size - 1) / bin_size;
}

inline HiCInteractionToBlockMapper::BlockMapperInter::BlockMapperInter(
    std::uint64_t block_bin_count, std::uint64_t block_column_count)
    : _block_bin_count(block_bin_count), _block_column_count(block_column_count) {
  assert(_block_bin_count != 0);
  assert(_block_column_count != 0);
}

inline std::uint64_t HiCInteractionToBlockMapper::BlockMapperInter::block_column_count() const {
  return _block_column_count;
}

inline std::uint64_t HiCInteractionToBlockMapper::BlockMapperInter::block_bin_count() const {
  return _block_bin_count;
}

inline std::uint64_t HiCInteractionToBlockMapper::BlockMapperInter::operator()(
    std::uint64_t bin1_id, std::uint64_t bin2_id) const {
  const auto i = bin1_id / block_bin_count();
  const auto j = bin2_id / block_bin_count();

  return (block_column_count() * j) + i;
}

inline HiCInteractionToBlockMapper::BlockMapperIntra::BlockMapperIntra(
    std::uint64_t block_bin_count, std::uint64_t block_column_count, std::int64_t base_depth)
    : _inter_mapper(block_bin_count, block_column_count), _base(init_base(base_depth)) {}

inline std::uint64_t HiCInteractionToBlockMapper::BlockMapperIntra::block_column_count() const {
  return _inter_mapper.block_column_count();
}

inline std::uint64_t HiCInteractionToBlockMapper::BlockMapperIntra::block_bin_count() const {
  return _inter_mapper.block_bin_count();
}

inline bool HiCInteractionToBlockMapper::BlockMapperIntra::use_inter_mapper() const noexcept {
  return _base == 0;
}

inline std::uint64_t HiCInteractionToBlockMapper::BlockMapperIntra::operator()(
    std::uint64_t bin1_id, std::uint64_t bin2_id) const {
  if (use_inter_mapper()) {
    return _inter_mapper(bin1_id, bin2_id);
  }
  const auto delta = bin1_id > bin2_id ? bin1_id - bin2_id : bin2_id - bin1_id;
  const auto n =
      static_cast<double>(delta) / std::sqrt(2.0) / static_cast<double>(block_bin_count());

  const auto depth = static_cast<std::uint64_t>(std::log(1.0 + n) / _base);
  const auto position_along_diagonal = (bin1_id + bin2_id) / 2 / block_bin_count();

  return depth * block_column_count() + position_along_diagonal;
}

inline double HiCInteractionToBlockMapper::BlockMapperIntra::init_base(
    std::int64_t base_depth) noexcept {
  if (base_depth > 1) {
    return std::log(static_cast<double>(base_depth));
  }
  if (base_depth < 0) {
    return static_cast<double>(-base_depth);
  }
  return std::log(2.0);
}

}  // namespace hictk::hic::internal
