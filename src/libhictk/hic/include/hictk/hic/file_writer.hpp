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

#include <BS_thread_pool.hpp>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <queue>
#include <string>

#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/binary_buffer.hpp"
#include "hictk/default_delete.hpp"
#include "hictk/expected_values_aggregator.hpp"
#include "hictk/filestream.hpp"
#include "hictk/hash.hpp"
#include "hictk/hic/file_writer_data_structures.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/hic/header.hpp"
#include "hictk/hic/interaction_block.hpp"
#include "hictk/hic/interaction_to_block_mapper.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk::hic::internal {

class HiCSectionOffsets {
  std::streampos _position{};
  std::size_t _size{};

 public:
  HiCSectionOffsets() = default;
  template <typename I1, typename I2>
  HiCSectionOffsets(I1 start_, I2 size_);

  [[nodiscard]] std::streampos start() const noexcept;
  [[nodiscard]] std::streampos end() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;

  void extend(std::size_t s) noexcept;
  void extend(std::streamoff s) noexcept;
  void set_size(std::size_t new_size) noexcept;
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
  struct Stats {
    double sum{};
    std::uint64_t nnz{};
  };

  filestream::FileStream<> _fs{};
  std::filesystem::path _tmpdir{};

  using BinTables = phmap::flat_hash_map<std::uint32_t, std::shared_ptr<const BinTable>>;
  using BlockIndex = phmap::btree_map<BlockIndexKey, phmap::btree_set<MatrixBlockMetadata>>;
  using BlockMappers = phmap::flat_hash_map<std::uint32_t, HiCInteractionToBlockMapper>;

  HiCHeader _header{};
  BinTables _bin_tables{};
  BlockIndex _block_index{};
  BlockMappers _block_mappers{};

  using StatsTank = phmap::flat_hash_map<std::uint32_t, Stats>;
  using FooterTank = phmap::btree_map<std::pair<Chromosome, Chromosome>, FooterMasterIndex>;

  MatrixBodyMetadataTank _matrix_metadata{};
  FooterTank _footers{};
  StatsTank _stats{};

  std::uint32_t _compression_lvl{};
  BinaryBuffer _bbuffer{};
  std::unique_ptr<libdeflate_compressor> _compressor{};
  std::string _compression_buffer{};

  phmap::btree_set<NormalizedExpectedValuesBlock> _normalized_expected_values{};
  phmap::btree_map<NormalizationVectorIndexBlock, std::vector<float>> _normalization_vectors{};

  HiCSectionOffsets _header_section{};
  HiCSectionOffsets _data_block_section{};
  HiCSectionOffsets _body_metadata_section{};
  HiCSectionOffsets _footer_section{};
  HiCSectionOffsets _expected_values_section{};
  HiCSectionOffsets _expected_values_norm_section{};
  HiCSectionOffsets _norm_vector_index_section{};
  HiCSectionOffsets _norm_vectors_section{};

  BS::thread_pool _tpool{};

  bool _skip_all_vs_all_matrix{};

  static constexpr std::uint32_t DEFAULT_CHROM_ALL_SCALE_FACTOR{1000};

 public:
  HiCFileWriter() = default;
  explicit HiCFileWriter(std::string_view path_, std::size_t n_threads = 1);
  HiCFileWriter(
      std::string_view path_, Reference chromosomes_, std::vector<std::uint32_t> resolutions_,
      std::string_view assembly_ = "unknown", std::size_t n_threads = 1,
      std::size_t chunk_size = 10'000'000,
      const std::filesystem::path& tmpdir = hictk::internal::TmpDir::default_temp_directory_path(),
      std::uint32_t compression_lvl = 11, bool skip_all_vs_all_matrix = false,
      std::size_t buffer_size = 32'000'000);

  [[nodiscard]] std::string_view path() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins(std::uint32_t resolution) const;
  [[nodiscard]] const std::vector<std::uint32_t>& resolutions() const noexcept;
  [[nodiscard]] auto stats(std::uint32_t resolution) const noexcept -> Stats;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void add_pixels(std::uint32_t resolution, PixelIt first_pixel, PixelIt last_pixel);

  // Write normalization vectors
  void add_norm_vector(std::string_view type, const Chromosome& chrom, std::string_view unit,
                       std::uint32_t bin_size, const balancing::Weights& weights,
                       bool force_overwrite = false,
                       std::size_t position = std::numeric_limits<std::size_t>::max(),
                       std::size_t n_bytes = std::numeric_limits<std::size_t>::max());
  void add_norm_vector(std::string_view type, std::string_view unit, std::uint32_t bin_size,
                       const balancing::Weights& weights, bool force_overwrite = false);

  void write_norm_vectors_and_norm_expected_values();

  void serialize();

 private:
  [[nodiscard]] static HiCHeader read_header(filestream::FileStream<>& fs);
  [[nodiscard]] static HiCHeader init_header(std::string_view path, Reference chromosomes,
                                             std::vector<std::uint32_t> resolutions,
                                             std::string_view assembly,
                                             bool skip_all_vs_all_matrix);
  [[nodiscard]] static auto init_bin_tables(const Reference& chromosomes,
                                            const std::vector<std::uint32_t>& resolutions)
      -> BinTables;
  [[nodiscard]] static auto init_interaction_block_mappers(const std::filesystem::path& root_folder,
                                                           const BinTables& bin_tables,
                                                           std::size_t chunk_size,
                                                           int compression_lvl) -> BlockMappers;
  [[nodiscard]] static BS::thread_pool init_tpool(std::size_t n_threads);

  // Write header
  void write_header();
  void write_footer_offset();
  void write_norm_vector_index();

  // Write pixels
  void write_pixels(bool skip_all_vs_all_matrix);
  HiCSectionOffsets write_pixels(const Chromosome& chrom1, const Chromosome& chrom2);
  HiCSectionOffsets write_pixels(const Chromosome& chrom1, const Chromosome& chrom2,
                                 std::uint32_t resolution);
  void write_all_matrix(std::uint32_t target_num_bins = 500);

  [[nodiscard]] HiCSectionOffsets write_interaction_block(
      std::streampos offset, std::uint64_t block_id, const Chromosome& chrom1,
      const Chromosome& chrom2, std::uint32_t resolution, const MatrixInteractionBlock<float>& blk);
  [[nodiscard]] auto write_interaction_blocks(std::streampos offset, const Chromosome& chrom1,
                                              const Chromosome& chrom2, std::uint32_t resolution)
      -> std::pair<HiCSectionOffsets, Stats>;

  // Normalization
  void add_norm_vector(const NormalizationVectorIndexBlock& blk, const balancing::Weights& weights,
                       bool force_overwrite = false);
  void add_norm_vector(const NormalizationVectorIndexBlock& blk, const std::vector<float>& weights,
                       bool force_overwrite = false);

  // Write body
  void write_body_metadata();
  void add_body_metadata(std::uint32_t resolution, const Chromosome& chrom1,
                         const Chromosome& chrom2, const std::string& unit = "BP");

  // Write footer
  void write_footers();
  void add_footer(const Chromosome& chrom1, const Chromosome& chrom2);
  void write_footer_size();

  HiCSectionOffsets write_empty_expected_values();
  HiCSectionOffsets write_empty_normalized_expected_values();
  HiCSectionOffsets compute_and_write_expected_values();
  HiCSectionOffsets compute_and_write_normalized_expected_values();
  HiCSectionOffsets write_norm_vectors();

  void finalize(bool compute_expected_values = false);

  [[nodiscard]] std::size_t compute_block_column_count(const Chromosome& chrom1,
                                                       const Chromosome& chrom2,
                                                       std::uint32_t resolution);
  [[nodiscard]] std::size_t compute_num_bins(const Chromosome& chrom1, const Chromosome& chrom2,
                                             std::uint32_t resolution);

  [[nodiscard]] ExpectedValuesBlock compute_expected_values(std::uint32_t resolution);
  [[nodiscard]] NormalizedExpectedValuesBlock compute_normalized_expected_values(
      std::uint32_t resolution, const balancing::Method& norm);

  void add_norm_expected_values(const NormalizedExpectedValuesBlock& blk,
                                bool force_overwrite = false);
  void read_norm_expected_values();
  void read_norm_vectors();
  [[nodiscard]] std::vector<float> read_norm_vector(const NormalizationVectorIndexBlock& blk);

  void read_offsets();

  // Methods to be called from worker threads
  auto merge_and_compress_blocks_thr(
      HiCInteractionToBlockMapper& mapper, std::mutex& mapper_mtx,
      std::queue<std::uint64_t>& block_id_queue, std::mutex& block_id_queue_mtx,
      moodycamel::BlockingConcurrentQueue<HiCInteractionToBlockMapper::BlockID>& block_queue,
      phmap::flat_hash_map<std::uint64_t, std::string>& serialized_block_tank,
      std::mutex& serialized_block_tank_mtx, std::atomic<bool>& early_return,
      std::uint64_t stop_token) -> Stats;
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
