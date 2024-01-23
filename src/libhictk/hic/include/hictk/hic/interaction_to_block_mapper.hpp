// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>
#if __has_include(<readerwriterqueue.h>)
#include <readerwriterqueue.h>
#else
#include <readerwriterqueue/readerwriterqueue.h>
#endif
#include <zstd.h>

#include <BS_thread_pool.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/default_delete.hpp"
#include "hictk/hic/binary_buffer.hpp"
#include "hictk/hic/file_writer_data_structures.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk::hic::internal {

template <typename N = float>
struct MatrixInteractionBlockFlat {
  std::vector<std::uint64_t> bin1_ids{};
  std::vector<std::uint64_t> bin2_ids{};
  std::vector<N> counts{};

  void emplace_back(ThinPixel<N>&& p);
  void emplace_back(Pixel<N>&& p);

  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, ZSTD_CCtx_s& compressor,
                                      std::string& compression_buffer, int compression_lvl,
                                      bool clear = true) const;
  [[nodiscard]] static std::vector<ThinPixel<N>> deserialize(BinaryBuffer& buffer,
                                                             ZSTD_DCtx_s& decompressor,
                                                             std::string& decompression_buffer);
};

class HiCInteractionToBlockMapper {
 public:
  class BlockMapperIntra;
  class BlockMapperInter;

  static constexpr std::uint32_t DEFAULT_INTRA_CUTOFF = 500;
  static constexpr std::uint32_t DEFAULT_INTER_CUTOFF = 5'000;
  static constexpr std::size_t DEFAULT_BLOCK_CAPACITY = 1'000;

  struct BlockID {
    std::uint32_t chrom1_id;
    std::uint32_t chrom2_id;
    std::uint64_t bid;

    [[nodiscard]] bool operator<(const BlockID& other) const noexcept;
    [[nodiscard]] bool operator==(const BlockID& other) const noexcept;
  };

  struct BlockIndex {
    std::uint64_t offset;
    std::uint32_t size;
  };

 private:
  std::filesystem::path _path{};
  filestream::FileStream _fs{};
  std::shared_ptr<const BinTable> _bin_table{};

  using BlockIndexMap = phmap::btree_map<BlockID, std::vector<BlockIndex>>;
  using ChromosomeIndexMap =
      phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, phmap::btree_set<BlockID>>;
  BlockIndexMap _block_index{};
  ChromosomeIndexMap _chromosome_index{};

  phmap::btree_map<BlockID, MatrixInteractionBlockFlat<float>> _blocks{};
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, float> _pixel_sums{};
  std::size_t _processed_pixels{};
  std::size_t _pending_pixels{};
  std::size_t _chunk_size{};

  phmap::flat_hash_map<Chromosome, BlockMapperIntra> _mappers_intra{};
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, BlockMapperInter> _mappers_inter{};

  BinaryBuffer _bbuffer{};
  int _compression_lvl{};
  std::unique_ptr<ZSTD_CCtx_s> _zstd_cctx{};
  std::unique_ptr<ZSTD_DCtx_s> _zstd_dctx{};
  std::string _compression_buffer{};

 public:
  HiCInteractionToBlockMapper() = default;
  HiCInteractionToBlockMapper(std::filesystem::path path, std::shared_ptr<const BinTable> bins,
                              std::size_t chunk_size, int compression_lvl);

  HiCInteractionToBlockMapper(const HiCInteractionToBlockMapper& other) = delete;
  HiCInteractionToBlockMapper(HiCInteractionToBlockMapper&& other) noexcept(noexcept_move_ctor()) =
      default;

  ~HiCInteractionToBlockMapper() noexcept;

  HiCInteractionToBlockMapper& operator=(const HiCInteractionToBlockMapper& other) = delete;
  HiCInteractionToBlockMapper& operator=(HiCInteractionToBlockMapper&& other) noexcept(
      noexcept_move_assignment_op()) = default;

  const Reference& chromosomes() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] bool empty(const Chromosome& chrom1, const Chromosome& chrom2) const noexcept;
  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(PixelIt first_pixel, PixelIt last_pixel,
                     std::uint32_t update_frequency = 10'000'000);
  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(PixelIt first_pixel, PixelIt last_pixel, BS::thread_pool& tpool,
                     std::uint32_t update_frequency = 10'000'000);

  [[nodiscard]] auto block_index() const noexcept -> const BlockIndexMap&;
  [[nodiscard]] auto chromosome_index() const noexcept -> const ChromosomeIndexMap&;
  [[nodiscard]] auto merge_blocks(const BlockID& bid) -> MatrixInteractionBlock<float>;
  [[nodiscard]] auto merge_blocks(const BlockID& bid, BinaryBuffer& bbuffer, ZSTD_DCtx_s& zstd_dctx,
                                  std::string& compression_buffer, std::mutex& mtx)
      -> MatrixInteractionBlock<float>;
  [[nodiscard]] float pixel_sum(const Chromosome& chrom1, const Chromosome& chrom2) const;
  [[nodiscard]] float pixel_sum() const;

  void finalize();
  void clear();

  [[nodiscard]] static std::size_t compute_block_column_count(
      const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t bin_size,
      std::uint32_t cutoff, std::size_t block_capacity = DEFAULT_BLOCK_CAPACITY);
  [[nodiscard]] static std::size_t compute_num_bins(const Chromosome& chrom1,
                                                    const Chromosome& chrom2, std::size_t bin_size);

 private:
  void init_block_mappers();

  template <typename N>
  [[nodiscard]] auto map(const ThinPixel<N>& p) const -> BlockID;
  template <typename N>
  [[nodiscard]] auto map(const Pixel<N>& p) const -> BlockID;

  template <typename N>
  void add_pixel(const ThinPixel<N>& p);
  template <typename N>
  void add_pixel(const Pixel<N>& p);

  [[nodiscard]] std::vector<Pixel<float>> fetch_pixels(const BlockID& bid);
  [[nodiscard]] std::vector<Pixel<float>> fetch_pixels(const BlockID& bid, BinaryBuffer& bbuffer,
                                                       ZSTD_DCtx_s& zstd_dctx,
                                                       std::string& compression_buffer,
                                                       std::mutex& mtx);

  void write_blocks();
  std::pair<std::uint64_t, std::uint32_t> write_block(const MatrixInteractionBlockFlat<float>& blk);

 public:
  class BlockMapperInter {
    std::uint64_t _block_bin_count{};
    std::uint64_t _block_column_count{};

   public:
    BlockMapperInter(std::uint64_t block_bin_count, std::uint64_t block_column_count);
    [[nodiscard]] std::uint64_t operator()(std::uint64_t bin1_id, std::uint64_t bin2_id) const;

    [[nodiscard]] std::uint64_t block_bin_count() const;
    [[nodiscard]] std::uint64_t block_column_count() const;
  };

  class BlockMapperIntra {
    BlockMapperInter _inter_mapper;
    double _base{};

    static constexpr std::int64_t DEFAULT_BASE_DEPTH = 2;

   public:
    BlockMapperIntra(std::uint64_t block_bin_count, std::uint64_t block_column_count,
                     std::int64_t base_depth = DEFAULT_BASE_DEPTH);
    [[nodiscard]] std::uint64_t operator()(std::uint64_t bin1_id, std::uint64_t bin2_id) const;

    [[nodiscard]] std::uint64_t block_bin_count() const;
    [[nodiscard]] std::uint64_t block_column_count() const;

   private:
    [[nodiscard]] bool use_inter_mapper() const noexcept;
    [[nodiscard]] static double init_base(std::int64_t base_depth) noexcept;
  };
};

}  // namespace hictk::hic::internal

template <>
struct std::hash<hictk::hic::internal::HiCInteractionToBlockMapper::BlockID> {
  inline std::size_t operator()(
      hictk::hic::internal::HiCInteractionToBlockMapper::BlockID const& bid) const noexcept {
    return hictk::internal::hash_combine(0, bid.chrom1_id, bid.chrom2_id, bid.bid);
  }
};

#include "./impl/interaction_to_block_mapper_impl.hpp"  // NOLINT
