// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>
#include <zstd.h>

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/binary_buffer.hpp"
#include "hictk/hic/file_writer_data_structures.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

template <>
struct std::default_delete<ZSTD_CCtx_s> {
  void operator()(ZSTD_CCtx_s* ctx) const { ZSTD_freeCCtx(ctx); }  // NOLINT
};

template <>
struct std::default_delete<ZSTD_DCtx_s> {
  void operator()(ZSTD_DCtx_s* ctx) const { ZSTD_freeDCtx(ctx); }  // NOLINT
};

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

  static constexpr std::int32_t DEFAULT_INTRA_CUTOFF = 500;
  static constexpr std::int32_t DEFAULT_INTER_CUTOFF = 5'000;
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
      phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, phmap::flat_hash_set<BlockID>>;
  BlockIndexMap _block_index{};
  ChromosomeIndexMap _chromosome_index{};

  phmap::btree_map<BlockID, MatrixInteractionBlockFlat<float>> _blocks{};
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, float> _pixel_sums{};
  std::size_t _processed_pixels{};
  std::size_t _pending_pixels{};

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
                              int compression_lvl);

  const Reference& chromosomes() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(PixelIt first_pixel, PixelIt last_pixel, std::size_t chunk_size = 100'000'000);

  [[nodiscard]] auto block_index() const noexcept -> const BlockIndexMap&;
  [[nodiscard]] auto chromosome_index() const noexcept -> const ChromosomeIndexMap&;
  [[nodiscard]] auto merge_blocks(const BlockID& bid) -> MatrixInteractionBlock<float>;
  [[nodiscard]] float pixel_sum(const Chromosome& chrom1, const Chromosome& chrom2) const;
  [[nodiscard]] float pixel_sum() const;

  void finalize();

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

  [[nodiscard]] std::vector<Pixel<float>> fetch_pixels(const BlockID& bid);

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
