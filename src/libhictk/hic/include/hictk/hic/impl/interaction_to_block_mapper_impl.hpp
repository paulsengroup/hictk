// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>
#include <zstd.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/binary_buffer.hpp"
#include "hictk/filestream.hpp"
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

inline bool HiCInteractionToBlockMapper::BlockID::operator>(const BlockID &other) const noexcept {
  if (chrom1_id != other.chrom1_id) {
    return chrom1_id > other.chrom1_id;
  }
  if (chrom2_id != other.chrom2_id) {
    return chrom2_id > other.chrom2_id;
  }
  return bid > other.bid;
}

inline bool HiCInteractionToBlockMapper::BlockID::operator==(const BlockID &other) const noexcept {
  return chrom1_id == other.chrom1_id && chrom2_id == other.chrom2_id && bid == other.bid;
}

inline HiCInteractionToBlockMapper::HiCInteractionToBlockMapper(
    std::filesystem::path path, std::shared_ptr<const BinTable> bins, std::size_t chunk_size,
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

inline HiCInteractionToBlockMapper::~HiCInteractionToBlockMapper() noexcept {
  try {
    _fs.close();
    std::filesystem::remove(_path);
  } catch (...) {
  }
}

inline const Reference &HiCInteractionToBlockMapper::chromosomes() const noexcept {
  return _bin_table->chromosomes();
}

inline std::size_t HiCInteractionToBlockMapper::size() const noexcept { return _processed_pixels; }
inline bool HiCInteractionToBlockMapper::empty() const noexcept { return size() == 0; }
inline bool HiCInteractionToBlockMapper::empty(const Chromosome &chrom1,
                                               const Chromosome &chrom2) const noexcept {
  auto it = _chromosome_index.find(std::make_pair(chrom1, chrom2));
  return it == _chromosome_index.end();
}

template <typename PixelIt, typename>
inline void HiCInteractionToBlockMapper::append_pixels(PixelIt first_pixel, PixelIt last_pixel,
                                                       std::uint32_t update_frequency) {
  using PixelT = remove_cvref_t<decltype(*first_pixel)>;
  static_assert(std::is_same_v<PixelT, ThinPixel<float>> || std::is_same_v<PixelT, Pixel<float>>);

  SPDLOG_DEBUG(FMT_STRING("mapping pixels to interaction blocks at resolution {}..."),
               _bin_table->resolution());

  auto t0 = std::chrono::steady_clock::now();
  for (std::size_t i = 0; first_pixel != last_pixel; ++i) {
    if (_pending_pixels >= _chunk_size) {
      write_blocks();
    }

    add_pixel(*first_pixel);
    std::ignore = ++first_pixel;

    if (i == update_frequency) {
      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      SPDLOG_INFO(FMT_STRING("ingesting pixels at {:.0f} pixels/s..."),
                  double(update_frequency) / delta);
      t0 = t1;
      i = 0;
    }
  }
}

namespace internal {
// I am declaring this as a freestanding funciton instead that as a lambda to work around a compiler
// bug in MSVC
template <typename PixelT>
Pixel<float> process_pixel_interaction_block(const BinTable &bin_table, PixelT &&pixel) {
  using PixelT_ = remove_cvref_t<PixelT>;
  static_assert(std::is_same_v<PixelT_, ThinPixel<float>> || std::is_same_v<PixelT_, Pixel<float>>);
  constexpr bool is_thin_pixel = std::is_same_v<PixelT_, ThinPixel<float>>;

  if constexpr (is_thin_pixel) {
    return Pixel(bin_table, pixel);
  } else {
    return pixel;
  }
}
}  // namespace internal

template <typename PixelIt, typename>
inline void HiCInteractionToBlockMapper::append_pixels(PixelIt first_pixel, PixelIt last_pixel,
                                                       BS::thread_pool &tpool,
                                                       std::uint32_t update_frequency) {
  if (tpool.get_thread_count() < 2) {
    return append_pixels(first_pixel, last_pixel);
  }

  SPDLOG_DEBUG(FMT_STRING("mapping pixels to interaction blocks using 2 threads..."));

  std::atomic<bool> early_return = false;
  moodycamel::BlockingReaderWriterQueue<Pixel<float>> queue(10'000);

  auto writer = tpool.submit_task([&]() {
    try {
      auto t0 = std::chrono::steady_clock::now();
      for (std::size_t i = 0; first_pixel != last_pixel && !early_return; ++i) {
        const auto pixel = internal::process_pixel_interaction_block(*_bin_table, *first_pixel);
        std::ignore = ++first_pixel;

        while (!queue.try_enqueue(pixel)) {
          if (early_return) {
            return;
          }
        }
        if (i == update_frequency) {
          const auto t1 = std::chrono::steady_clock::now();
          const auto delta =
              static_cast<double>(
                  std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
              1000.0;
          SPDLOG_INFO(FMT_STRING("ingesting pixels at {:.0f} pixels/s..."),
                      double(update_frequency) / delta);
          t0 = t1;
          i = 0;
        }
      }
      queue.enqueue(Pixel<float>{});

    } catch (...) {
      early_return = true;
      throw;
    }
  });

  auto reader = tpool.submit_task([&]() {
    try {
      Pixel<float> p{};
      while (!early_return) {
        if (_pending_pixels >= _chunk_size) {
          write_blocks();
        }

        queue.wait_dequeue(p);
        if (p.count == 0) {
          return;
        }

        add_pixel(p);
        ++_processed_pixels;
        ++_pending_pixels;
      }
    } catch (...) {
      early_return = true;
      throw;
    }
  });

  writer.get();
  reader.get();
}

inline auto HiCInteractionToBlockMapper::block_index() const noexcept -> const BlockIndexMap & {
  return _block_index;
}

inline auto HiCInteractionToBlockMapper::chromosome_index() const noexcept
    -> const MatrixIndexMap & {
  return _chromosome_index;
}

inline auto HiCInteractionToBlockMapper::merge_blocks(const BlockID &bid)
    -> MatrixInteractionBlock<float> {
  std::mutex dummy_mtx{};
  return merge_blocks(bid, _bbuffer, *_zstd_dctx, _compression_buffer, dummy_mtx);
}

inline auto HiCInteractionToBlockMapper::merge_blocks(const BlockID &bid, BinaryBuffer &bbuffer,
                                                      ZSTD_DCtx_s &zstd_dctx,
                                                      std::string &compression_buffer,
                                                      std::mutex &mtx)
    -> MatrixInteractionBlock<float> {
  MatrixInteractionBlock<float> blk{};
  for (auto &&pixel : fetch_pixels(bid, bbuffer, zstd_dctx, compression_buffer, mtx)) {
    blk.emplace_back(std::move(pixel));
  }
  blk.finalize();
  return blk;
}

inline float HiCInteractionToBlockMapper::pixel_sum(const Chromosome &chrom1,
                                                    const Chromosome &chrom2) const {
  auto match = _pixel_sums.find(std::make_pair(chrom1, chrom2));
  if (match == _pixel_sums.end()) {
    return 0;
  }
  return match->second;
}

inline float HiCInteractionToBlockMapper::pixel_sum() const {
  return std::accumulate(
      _pixel_sums.begin(), _pixel_sums.end(), 0.0F,
      [](const float accumulator, const auto &kv) { return accumulator + kv.second; });
}

inline void HiCInteractionToBlockMapper::finalize() {
  if (_processed_pixels > _pending_pixels) {
    write_blocks();
  }
}

inline void HiCInteractionToBlockMapper::clear() {
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

inline void HiCInteractionToBlockMapper::init_block_mappers() {
  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);

      const auto num_bins = compute_num_bins(chrom1, chrom2, _bin_table->resolution());
      const auto num_columns = compute_block_column_count(
          chrom1, chrom2, _bin_table->resolution(),
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

template <typename N>
inline void HiCInteractionToBlockMapper::add_pixel(const ThinPixel<N> &p) {
  add_pixel(Pixel(*_bin_table, p));
}

template <typename N>
inline void HiCInteractionToBlockMapper::add_pixel(const Pixel<N> &p) {
  auto bid = map(p);

  const auto &chrom1 = p.coords.bin1.chrom();
  const auto &chrom2 = p.coords.bin2.chrom();
  const auto chrom_pair = std::make_pair(chrom1, chrom2);

  auto match1 = _blocks.find(bid);
  if (match1 != _blocks.end()) {
    _pixel_sums.at(chrom_pair) += p.count;
    match1->second.emplace_back(p.to_thin());
  } else {
    _pixel_sums.emplace(chrom_pair, p.count);
    auto [it, _] = _blocks.emplace(std::move(bid), MatrixInteractionBlockFlat<float>{});
    it->second.emplace_back(p.to_thin());
  }

  auto match2 = _chromosome_index.find(chrom_pair);
  if (match2 != _chromosome_index.end()) {
    match2->second.emplace(bid);
  } else {
    _chromosome_index.emplace(chrom_pair, phmap::btree_set<BlockID>{bid});
  }
  ++_processed_pixels;
  ++_pending_pixels;
}

inline std::vector<Pixel<float>> HiCInteractionToBlockMapper::fetch_pixels(const BlockID &bid) {
  std::mutex dummy_mtx{};
  return fetch_pixels(bid, _bbuffer, *_zstd_dctx, _compression_buffer, dummy_mtx);
}

inline std::vector<Pixel<float>> HiCInteractionToBlockMapper::fetch_pixels(
    const BlockID &bid, BinaryBuffer &bbuffer, ZSTD_DCtx_s &zstd_dctx,
    std::string &compression_buffer, std::mutex &mtx) {
  std::vector<Pixel<float>> pixels{};
  auto match = _blocks.find(bid);
  if (match != _blocks.end()) {
    const auto &flat_pixels = match->second;
    pixels.reserve(flat_pixels.size());
    for (std::size_t i = 0; i < flat_pixels.size(); ++i) {
      pixels.emplace_back(
          Pixel(*_bin_table, ThinPixel<float>{flat_pixels.bin1_ids[i], flat_pixels.bin2_ids[i],
                                              flat_pixels.counts[i]}));
    }
    return pixels;
  }

  for (const auto &[pos, size] : _block_index.at(bid)) {
    {
      std::scoped_lock lck(mtx);
      _fs.seekg(static_cast<std::streamoff>(pos));
      _fs.read(bbuffer.reset(), size);
    }
    const auto flat_pixels =
        MatrixInteractionBlockFlat<float>::deserialize(bbuffer, zstd_dctx, compression_buffer);
    pixels.reserve(pixels.size() + flat_pixels.size());
    for (const auto &p : flat_pixels) {
      pixels.emplace_back(Pixel(*_bin_table, p));
    }
  }
  return pixels;
}

inline void HiCInteractionToBlockMapper::write_blocks() {
  if (!std::filesystem::exists(_path)) {
    _fs = filestream::FileStream<>::create(_path.string(), std::make_shared<std::mutex>());
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
