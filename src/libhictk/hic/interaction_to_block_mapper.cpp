// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/interaction_to_block_mapper.hpp"

#include <fmt/std.h>
#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>
#include <zstd.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iosfwd>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/binary_buffer.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/filestream.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

bool HiCInteractionToBlockMapper::BlockID::operator<(const BlockID &other) const noexcept {
  if (chrom1_id != other.chrom1_id) {
    return chrom1_id < other.chrom1_id;
  }
  if (chrom2_id != other.chrom2_id) {
    return chrom2_id < other.chrom2_id;
  }
  return bid < other.bid;
}

bool HiCInteractionToBlockMapper::BlockID::operator>(const BlockID &other) const noexcept {
  if (chrom1_id != other.chrom1_id) {
    return chrom1_id > other.chrom1_id;
  }
  if (chrom2_id != other.chrom2_id) {
    return chrom2_id > other.chrom2_id;
  }
  return bid > other.bid;
}

bool HiCInteractionToBlockMapper::BlockID::operator==(const BlockID &other) const noexcept {
  return chrom1_id == other.chrom1_id && chrom2_id == other.chrom2_id && bid == other.bid;
}

HiCInteractionToBlockMapper::HiCInteractionToBlockMapper(std::filesystem::path path,
                                                         std::shared_ptr<const BinTable> bins,
                                                         std::size_t chunk_size,
                                                         int compression_lvl)
    : _path(std::move(path)),
      _bin_table(std::move(bins)),
      _chunk_size(chunk_size),
      _compression_lvl(compression_lvl),
      _zstd_cctx(ZSTD_createCCtx()),
      _zstd_dctx(ZSTD_createDCtx()) {
  assert(_chunk_size != 0);
  SPDLOG_INFO(
      FMT_STRING("initializing HiCInteractionToBlockMapper using \"{}\" as temporary file..."),
      _path);
  init_block_mappers();
}

// NOLINTNEXTLINE(bugprone-exception-escape)
HiCInteractionToBlockMapper::~HiCInteractionToBlockMapper() noexcept {
  try {
    _fs.close();
    std::filesystem::remove(_path);
  } catch (const std::exception &e) {
    SPDLOG_WARN(FMT_STRING("failed to remove temporary file \"{}\": {}"), _path, e.what());
  } catch (...) {
    SPDLOG_WARN(FMT_STRING("failed to remove temporary file \"{}\""), _path);
  }
}

const Reference &HiCInteractionToBlockMapper::chromosomes() const noexcept {
  return _bin_table->chromosomes();
}

std::size_t HiCInteractionToBlockMapper::size() const noexcept { return _processed_pixels; }
bool HiCInteractionToBlockMapper::empty() const noexcept { return size() == 0; }
bool HiCInteractionToBlockMapper::empty(const Chromosome &chrom1,
                                        const Chromosome &chrom2) const noexcept {
  auto it = _chromosome_index.find(std::make_pair(chrom1, chrom2));
  return it == _chromosome_index.end();
}

auto HiCInteractionToBlockMapper::block_index() const noexcept -> const BlockIndexMap & {
  return _block_index;
}

auto HiCInteractionToBlockMapper::chromosome_index() const noexcept -> const MatrixIndexMap & {
  return _chromosome_index;
}

auto HiCInteractionToBlockMapper::merge_blocks(const BlockID &bid)
    -> MatrixInteractionBlock<float> {
  std::mutex dummy_mtx{};
  return merge_blocks(bid, _bbuffer, *_zstd_dctx, _compression_buffer, dummy_mtx);
}

auto HiCInteractionToBlockMapper::merge_blocks(const BlockID &bid, BinaryBuffer &bbuffer,
                                               ZSTD_DCtx_s &zstd_dctx,
                                               std::string &compression_buffer, std::mutex &mtx)
    -> MatrixInteractionBlock<float> {
  MatrixInteractionBlock<float> blk{};
  for (const auto &pixel : fetch_pixels(bid, bbuffer, zstd_dctx, compression_buffer, mtx)) {
    blk.emplace_back(pixel);
  }
  blk.finalize();
  return blk;
}

float HiCInteractionToBlockMapper::pixel_sum(const Chromosome &chrom1,
                                             const Chromosome &chrom2) const {
  auto match = _pixel_sums.find(std::make_pair(chrom1, chrom2));
  if (match == _pixel_sums.end()) {
    return 0;
  }
  return match->second;
}

float HiCInteractionToBlockMapper::pixel_sum() const {
  return std::accumulate(
      _pixel_sums.begin(), _pixel_sums.end(), 0.0F,
      [](const float accumulator, const auto &kv) { return accumulator + kv.second; });
}

void HiCInteractionToBlockMapper::finalize() {
  if (_processed_pixels > _pending_pixels) {
    write_blocks();
  }
}

void HiCInteractionToBlockMapper::clear() {
  _blocks.clear();
  _block_index.clear();
  _chromosome_index.clear();
  _pixel_sums.clear();
  _processed_pixels = 0;
  _pending_pixels = 0;
  _bbuffer.reset().shrink_to_fit();
  _compression_buffer.clear();
  _compression_buffer.shrink_to_fit();
  std::error_code ec{};
  std::filesystem::remove(_path, ec);
}

void HiCInteractionToBlockMapper::init_block_mappers() {
  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);

      const auto num_bins = compute_num_bins(chrom1, chrom2, _bin_table->resolution());
      const auto num_columns = compute_block_column_count(
          chrom1, chrom2, _bin_table->resolution(),
          chrom1 == chrom2 ? DEFAULT_INTRA_CUTOFF : DEFAULT_INTER_CUTOFF);
      const auto num_rows = (num_bins / num_columns) + 1;

      if (chrom1 == chrom2) {
        _mappers_intra.emplace(chrom1, BlockMapperIntra{num_rows, num_columns});
      } else {
        _mappers_inter.emplace(std::make_pair(chrom1, chrom2),
                               BlockMapperInter{num_rows, num_columns});
      }
    }
  }
}

std::pair<std::uint64_t, std::uint32_t> HiCInteractionToBlockMapper::write_block(
    const MatrixInteractionBlockFlat<float> &blk) {
  const auto offset = _fs.tellp();
  _fs.write(blk.serialize(_bbuffer, *_zstd_cctx, _compression_buffer, _compression_lvl));
  const auto size = _fs.tellp() - offset;
  return std::make_pair(offset, static_cast<std::uint32_t>(size));
}

std::vector<Pixel<float>> HiCInteractionToBlockMapper::fetch_pixels(const BlockID &bid) {
  std::mutex dummy_mtx{};
  return fetch_pixels(bid, _bbuffer, *_zstd_dctx, _compression_buffer, dummy_mtx);
}

std::vector<Pixel<float>> HiCInteractionToBlockMapper::fetch_pixels(const BlockID &bid,
                                                                    BinaryBuffer &bbuffer,
                                                                    ZSTD_DCtx_s &zstd_dctx,
                                                                    std::string &compression_buffer,
                                                                    std::mutex &mtx) {
  std::vector<Pixel<float>> pixels{};
  auto match = _blocks.find(bid);
  if (match != _blocks.end()) {
    const auto &flat_pixels = match->second;
    pixels.reserve(flat_pixels.size());
    for (std::size_t i = 0; i < flat_pixels.size(); ++i) {
      pixels.emplace_back(*_bin_table,
                          ThinPixel<float>{flat_pixels.bin1_ids[i], flat_pixels.bin2_ids[i],
                                           flat_pixels.counts[i]});
    }
    return pixels;
  }

  for (const auto &[pos, size] : _block_index.at(bid)) {
    {
      [[maybe_unused]] const std::scoped_lock lck(mtx);
      _fs.seekg(static_cast<std::streamoff>(pos));
      _fs.read(bbuffer.reset(), size);
    }
    const auto flat_pixels =
        MatrixInteractionBlockFlat<float>::deserialize(bbuffer, zstd_dctx, compression_buffer);
    pixels.reserve(pixels.size() + flat_pixels.size());
    for (const auto &p : flat_pixels) {
      pixels.emplace_back(*_bin_table, p);
    }
  }
  return pixels;
}

void HiCInteractionToBlockMapper::write_blocks() {
  if (!std::filesystem::exists(_path)) {
    _fs = filestream::FileStream::create(_path.string(), std::make_shared<std::mutex>());
  }
  SPDLOG_DEBUG(FMT_STRING("writing {} pixels to file {}..."), _pending_pixels, _path);
  for (auto &[bid, blk] : _blocks) {
    const auto [offset, size] = write_block(blk);
    auto match = _block_index.find(bid);
    if (match != _block_index.end()) {
      match->second.emplace_back(BlockIndex{offset, size});
    } else {
      _block_index.emplace(bid, std::vector<BlockIndex>{{offset, size}});
    }
  }
  _fs.flush();
  _blocks.clear();
  _pending_pixels = 0;
}

std::size_t HiCInteractionToBlockMapper::compute_block_column_count(const Chromosome &chrom1,
                                                                    const Chromosome &chrom2,
                                                                    std::uint32_t bin_size,
                                                                    std::uint32_t cutoff,
                                                                    std::size_t block_capacity) {
  const auto num_bins = compute_num_bins(chrom1, chrom2, bin_size);
  auto num_columns = (num_bins / block_capacity) + 1;
  if (bin_size < cutoff) {
    const auto genome_size = num_bins * bin_size;
    num_columns = genome_size / (block_capacity * cutoff);
  }

  const auto max_sqrt =
      static_cast<std::size_t>(std::sqrt(std::numeric_limits<std::int32_t>::max()));
  return std::clamp(num_columns, std::size_t{1}, max_sqrt - 1);
}

std::size_t HiCInteractionToBlockMapper::compute_num_bins(const Chromosome &chrom1,
                                                          const Chromosome &chrom2,
                                                          std::size_t bin_size) {
  const auto max_size = std::max(chrom1.size(), chrom2.size());
  return (max_size + bin_size - 1) / bin_size;
}

HiCInteractionToBlockMapper::BlockMapperInter::BlockMapperInter(std::uint64_t block_bin_count,
                                                                std::uint64_t block_column_count)
    : _block_bin_count(block_bin_count), _block_column_count(block_column_count) {
  assert(_block_bin_count != 0);
  assert(_block_column_count != 0);
}

std::uint64_t HiCInteractionToBlockMapper::BlockMapperInter::block_column_count() const {
  return _block_column_count;
}

std::uint64_t HiCInteractionToBlockMapper::BlockMapperInter::block_bin_count() const {
  return _block_bin_count;
}

std::uint64_t HiCInteractionToBlockMapper::BlockMapperInter::operator()(
    std::uint64_t bin1_id, std::uint64_t bin2_id) const {
  const auto i = bin1_id / block_bin_count();
  const auto j = bin2_id / block_bin_count();

  return (block_column_count() * j) + i;
}

HiCInteractionToBlockMapper::BlockMapperIntra::BlockMapperIntra(std::uint64_t block_bin_count,
                                                                std::uint64_t block_column_count,
                                                                std::int64_t base_depth)
    : _inter_mapper(block_bin_count, block_column_count), _base(init_base(base_depth)) {}

std::uint64_t HiCInteractionToBlockMapper::BlockMapperIntra::block_column_count() const {
  return _inter_mapper.block_column_count();
}

std::uint64_t HiCInteractionToBlockMapper::BlockMapperIntra::block_bin_count() const {
  return _inter_mapper.block_bin_count();
}

bool HiCInteractionToBlockMapper::BlockMapperIntra::use_inter_mapper() const noexcept {
  return _base == 0;
}

std::uint64_t HiCInteractionToBlockMapper::BlockMapperIntra::operator()(
    std::uint64_t bin1_id, std::uint64_t bin2_id) const {
  if (use_inter_mapper()) {
    return _inter_mapper(bin1_id, bin2_id);
  }
  const auto delta = bin1_id > bin2_id ? bin1_id - bin2_id : bin2_id - bin1_id;
  const auto n =
      static_cast<double>(delta) / std::sqrt(2.0) / static_cast<double>(block_bin_count());

  const auto depth = static_cast<std::uint64_t>(std::log(1.0 + n) / _base);
  const auto position_along_diagonal = (bin1_id + bin2_id) / 2 / block_bin_count();

  return (depth * block_column_count()) + position_along_diagonal;
}

double HiCInteractionToBlockMapper::BlockMapperIntra::init_base(std::int64_t base_depth) noexcept {
  if (base_depth > 1) {
    return std::log(static_cast<double>(base_depth));
  }
  if (base_depth < 0) {
    return static_cast<double>(-base_depth);
  }
  return std::log(2.0);  // NOLINT(*-avoid-magic-numbers)
}

}  // namespace hictk::hic::internal
