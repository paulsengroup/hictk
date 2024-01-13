// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#if __has_include(<blockingconcurrentqueue.h>)
#include <blockingconcurrentqueue.h>
#else
#include <concurrentqueue/blockingconcurrentqueue.h>
#endif
#include <libdeflate.h>
#include <parallel_hashmap/btree.h>
#include <zstd.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <queue>
#include <string>

#include "hictk/bin_table.hpp"
#include "hictk/hash.hpp"
#include "hictk/hic/binary_buffer.hpp"
#include "hictk/hic/expected_values_aggregator.hpp"
#include "hictk/hic/file_writer_data_structures.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/hic/header.hpp"
#include "hictk/hic/interaction_block.hpp"
#include "hictk/hic/interaction_to_block_mapper.hpp"
#include "hictk/tmpdir.hpp"

template <>
struct std::default_delete<libdeflate_compressor> {
  void operator()(libdeflate_compressor* compressor) const {
    libdeflate_free_compressor(compressor);
  }
};

namespace hictk::hic::internal {

class HiCSectionOffsets {
  std::streamoff _position{};
  std::size_t _size{};

 public:
  HiCSectionOffsets() = default;
  template <typename I1, typename I2>
  HiCSectionOffsets(I1 start_, I2 size_);

  [[nodiscard]] std::streamoff start() const noexcept;
  [[nodiscard]] std::streamoff end() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::size_t& size() noexcept;
};

struct BlockIndexKey {
  Chromosome chrom1;
  Chromosome chrom2;
  std::uint32_t resolution;

  [[nodiscard]] bool operator<(const BlockIndexKey& other) const noexcept;
};

class MatrixBodyMetadataTank {
 public:
  struct Key {
    Chromosome chrom1{};
    Chromosome chrom2{};

    bool operator==(const Key& other) const noexcept;
  };

 private:
  phmap::flat_hash_map<Key, MatrixBodyMetadata> _tank{};
  phmap::flat_hash_map<Key, HiCSectionOffsets> _offsets{};

 public:
  MatrixBodyMetadataTank() = default;

  [[nodiscard]] bool contains(const Chromosome& chrom1, const Chromosome& chrom2) const noexcept;
  [[nodiscard]] auto at(const Chromosome& chrom1, const Chromosome& chrom2) const
      -> const MatrixBodyMetadata&;
  [[nodiscard]] HiCSectionOffsets offset(const Chromosome& chrom1, const Chromosome& chrom2) const;

  void insert(const Chromosome& chrom1, const Chromosome& chrom2, MatrixMetadata matrix_metadata,
              MatrixResolutionMetadata matrix_resolution_metadata);
  void update_offsets(const Chromosome& chrom1, const Chromosome& chrom2, std::streamoff position,
                      std::size_t size);

  void remove(const Chromosome& chrom1, const Chromosome& chrom2);

  auto operator()() const noexcept -> const phmap::flat_hash_map<Key, MatrixBodyMetadata>&;
};

class HiCFileWriter {
  std::shared_ptr<const HiCHeader> _header{};
  std::shared_ptr<filestream::FileStream> _fs{};
  std::unique_ptr<const hictk::internal::TmpDir> _tmpdir{};

  using BinTables = phmap::flat_hash_map<std::uint32_t, std::shared_ptr<const BinTable>>;
  using BlockIndex = phmap::btree_map<BlockIndexKey, phmap::btree_set<MatrixBlockMetadata>>;
  using BlockMappers = phmap::flat_hash_map<std::uint32_t, HiCInteractionToBlockMapper>;

  BinTables _bin_tables{};
  BlockIndex _block_index{};
  BlockMappers _block_mappers{};

  MatrixBodyMetadataTank _matrix_metadata{};
  using FooterTank = phmap::btree_map<std::pair<Chromosome, Chromosome>, FooterV5>;
  FooterTank _footers{};

  std::int32_t _compression_lvl{};
  BinaryBuffer _bbuffer{};
  std::unique_ptr<libdeflate_compressor> _compressor{};
  std::string _compression_buffer{};

  HiCSectionOffsets _header_section{};
  HiCSectionOffsets _data_block_section{};
  HiCSectionOffsets _body_metadata_section{};
  HiCSectionOffsets _footer_section{};
  HiCSectionOffsets _expected_values_section{};
  HiCSectionOffsets _expected_values_norm_section{};
  HiCSectionOffsets _norm_vectors_section{};

  BS::thread_pool _tpool{};

 public:
  HiCFileWriter() = default;
  explicit HiCFileWriter(
      HiCHeader header, std::size_t n_threads = 1, std::size_t chunk_size = 10'000'000,
      const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
      std::int32_t compression_lvl = 9, std::size_t buffer_size = 32'000'000);

  [[nodiscard]] std::string_view url() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins(std::uint32_t resolution) const;
  [[nodiscard]] const std::vector<std::uint32_t>& resolutions() const noexcept;

  void serialize();

  // Write header
  void write_header();
  void write_footer_offset();
  void write_norm_vector_index();

  // Write pixels
  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void add_pixels(PixelIt first_pixel, PixelIt last_pixel);
  void write_pixels();
  void write_all_matrix(std::uint32_t target_resolution = 2'500'000);

  // Write body
  void write_body_metadata();
  void add_body_metadata(std::uint32_t resolution, const Chromosome& chrom1,
                         const Chromosome& chrom2, const std::string& unit = "BP");

  // Write footer
  void write_footers();
  void add_footer(const Chromosome& chrom1, const Chromosome& chrom2);
  void write_footer_size();

  // Write expected/normalization values
  void compute_and_write_expected_values();

  void write_empty_expected_values();
  void write_empty_normalized_expected_values();
  void write_empty_norm_vectors();

  void finalize();

 private:
  [[nodiscard]] static std::shared_ptr<const HiCHeader> init_header(
      HiCHeader&& header, std::uint32_t all_scale_factor = 1);
  [[nodiscard]] static auto init_bin_tables(const Reference& chromosomes,
                                            const std::vector<std::uint32_t>& resolutions)
      -> BinTables;
  [[nodiscard]] static auto init_interaction_block_mappers(const std::filesystem::path& root_folder,
                                                           const BinTables& bin_tables,
                                                           std::size_t chunk_size,
                                                           int compression_lvl) -> BlockMappers;
  [[nodiscard]] BS::thread_pool init_tpool(std::size_t n_threads);

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void add_pixels(std::uint32_t resolution, PixelIt first_pixel, PixelIt last_pixel);

  auto write_pixels(const Chromosome& chrom1, const Chromosome& chrom2) -> HiCSectionOffsets;
  auto write_pixels(const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t resolution)
      -> HiCSectionOffsets;

  auto write_interaction_block(std::uint64_t block_id, const Chromosome& chrom1,
                               const Chromosome& chrom2, std::uint32_t resolution,
                               const MatrixInteractionBlock<float>& blk) -> HiCSectionOffsets;
  std::size_t write_interaction_blocks(const Chromosome& chrom1, const Chromosome& chrom2,
                                       std::uint32_t resolution);

  [[nodiscard]] std::size_t compute_block_column_count(const Chromosome& chrom1,
                                                       const Chromosome& chrom2,
                                                       std::uint32_t resolution);
  [[nodiscard]] std::size_t compute_num_bins(const Chromosome& chrom1, const Chromosome& chrom2,
                                             std::uint32_t resolution);

  // Methods to be called from worker threads
  std::size_t merge_and_compress_blocks_thr(
      HiCInteractionToBlockMapper& mapper, std::mutex& mapper_mtx,
      std::queue<std::uint64_t>& block_id_queue, std::mutex& block_id_queue_mtx,
      moodycamel::BlockingConcurrentQueue<HiCInteractionToBlockMapper::BlockID>& block_queue,
      phmap::flat_hash_map<std::uint64_t, std::string>& serialized_block_tank,
      std::mutex& serialized_block_tank_mtx, std::atomic<bool>& early_return,
      std::uint64_t stop_token);
  void write_compressed_blocks_thr(
      const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t resolution,
      std::queue<std::uint64_t>& block_id_queue, std::mutex& block_id_queue_mtx,
      phmap::flat_hash_map<std::uint64_t, std::string>& serialized_block_tank,
      std::mutex& serialized_block_tank_mtx, std::atomic<bool>& early_return,
      std::uint64_t stop_token);
};
}  // namespace hictk::hic::internal

template <>
struct std::hash<hictk::hic::internal::MatrixBodyMetadataTank::Key> {
  inline std::size_t operator()(
      hictk::hic::internal::MatrixBodyMetadataTank::Key const& k) const noexcept {
    return hictk::internal::hash_combine(0, k.chrom1, k.chrom2);
  }
};

#include "./impl/file_writer_impl.hpp"  // NOLINT
