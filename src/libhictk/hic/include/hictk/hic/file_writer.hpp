// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <libdeflate.h>
#include <parallel_hashmap/btree.h>
#include <zstd.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

#include "hictk/bin_table.hpp"
#include "hictk/hic/binary_buffer.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/hic/header.hpp"
#include "hictk/hic/interaction_block.hpp"
#include "hictk/tmpdir.hpp"

template <>
struct std::default_delete<libdeflate_compressor> {
  void operator()(libdeflate_compressor* compressor) const {
    libdeflate_free_compressor(compressor);
  }
};

template <>
struct std::default_delete<ZSTD_CCtx_s> {
  void operator()(ZSTD_CCtx_s* ctx) const { ZSTD_freeCCtx(ctx); }  // NOLINT
};

template <>
struct std::default_delete<ZSTD_DCtx_s> {
  void operator()(ZSTD_DCtx_s* ctx) const { ZSTD_freeDCtx(ctx); }  // NOLINT
};

namespace hictk::hic::internal {

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#matrix-metadata
struct MatrixMetadata {
  std::int32_t chr1Idx{};
  std::int32_t chr2Idx{};
  std::int32_t nResolutions{};

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

struct MatrixBlockMetadata {
  std::int32_t blockNumber{};
  std::int64_t blockPosition{};
  std::int32_t blockSizeBytes{};

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
  [[nodiscard]] bool operator<(const MatrixBlockMetadata& other) const noexcept;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#resolution-zoom-level-metadata
struct MatrixResolutionMetadata {
  std::string unit{};
  std::int32_t resIdx{};
  float sumCounts{};
  std::int32_t occupiedCellCount = 0;  // Not used
  float percent5 = 0;                  // Not used
  float percent95 = 0;                 // Not used
  std::int32_t binSize{};
  std::int32_t blockSize{};
  std::int32_t blockColumnCount{};
  std::int32_t blockCount{};

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#blocks
template <typename N = float>
struct MatrixInteractionBlock {
  std::int32_t nRecords{};
  std::int32_t binColumnOffset{std::numeric_limits<std::int32_t>::max()};
  std::int32_t binRowOffset{std::numeric_limits<std::int32_t>::max()};
  std::uint8_t useFloatContact{};
  std::uint8_t useIntXPos{};
  std::uint8_t useIntYPos{};
  std::uint8_t matrixRepresentation{};

  void emplace_back(Pixel<N>&& p);
  void finalize();

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, libdeflate_compressor& compressor,
                                      std::string& compression_buffer, bool clear = true) const;

 private:
  using RowID = std::int32_t;
  using Row = std::vector<Pixel<N>>;
  phmap::btree_map<RowID, Row> _interactions;
};

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

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#master-index
struct MasterIndex {
  std::string key;
  std::int64_t position;
  std::int32_t size;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

struct ExpectedValuesBlock {
  std::string unit{};
  std::int32_t binSize{};
  std::int64_t nValues{};
  std::vector<float> value{};
  std::int32_t nChrScaleFactors{};
  std::vector<std::int32_t> chrIndex{};
  std::vector<float> chrScaleFactor{};

  ExpectedValuesBlock(std::string_view unit_, std::uint32_t bin_size,
                      const std::vector<double>& weights,
                      const std::vector<std::uint32_t>& chrom_ids,
                      const std::vector<double>& scale_factors);
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#expected-value-vectors
struct ExpectedValues {
  std::int32_t nExpectedValueVectors = 0;
  std::vector<ExpectedValuesBlock> expectedValues;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalized-expected-value-vectors
struct NormalizedExpectedValues {
  std::int32_t nNormExpectedValueVectors = 0;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalization-vector-index
struct NormalizationVectorIndex {
  std::int32_t nNormVectors = 0;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalization-vector-arrays-1-per-normalization-vector
struct NormalizationVectorArray {
  std::int64_t nValues = 0;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

struct FooterV5 {
  MasterIndex masterIndex{};

  ExpectedValues expectedValues{};
  NormalizedExpectedValues normExpectedValues{};
  NormalizationVectorIndex normVectIndex{};
  std::vector<NormalizationVectorArray> normVectArray{};

  FooterV5() = default;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

struct BlockIndexKey {
  Chromosome chrom1;
  Chromosome chrom2;
  std::uint32_t resolution;

  [[nodiscard]] bool operator<(const BlockIndexKey& other) const noexcept;
};

class ExpectedValuesAggregator {
  std::shared_ptr<const BinTable> _bins{};
  std::size_t _num_bins_gw{};

  using CisKey = Chromosome;
  using TransKey = std::pair<CisKey, CisKey>;
  phmap::flat_hash_map<CisKey, double> _cis_sum{};
  phmap::flat_hash_map<TransKey, double> _trans_sum{};

  std::vector<double> _possible_distances{};
  std::vector<double> _actual_distances{};

  std::vector<double> _weights{};
  phmap::btree_map<Chromosome, double> _scaling_factors{};

 public:
  ExpectedValuesAggregator() = default;
  explicit ExpectedValuesAggregator(std::shared_ptr<const BinTable> bins);
  void add(const ThinPixel<float>& p);
  void add(const Pixel<float>& p);

  void compute_density();

  [[nodiscard]] const std::vector<double>& weights() const noexcept;

  [[nodiscard]] double scaling_factor(const Chromosome& chrom) const;
  [[nodiscard]] const phmap::btree_map<Chromosome, double>& scaling_factors() const noexcept;

 private:
  [[nodiscard]] const Reference& chromosomes() const noexcept;

  void compute_density_cis();
  void compute_density_trans();

  [[nodiscard]] double at(const Chromosome& chrom) const;
  [[nodiscard]] double at(const Chromosome& chrom1, const Chromosome& chrom2) const;

  [[nodiscard]] double& at(const Chromosome& chrom);
  [[nodiscard]] double& at(const Chromosome& chrom1, const Chromosome& chrom2);
};

template <typename N>
class PixelTank {
  using ChromosomePair = std::pair<Chromosome, Chromosome>;
  using ChromPixelTank = std::vector<ThinPixel<N>>;

  std::shared_ptr<const BinTable> _bin_table{};
  phmap::flat_hash_map<ChromosomePair, ChromPixelTank> _pixel_tank{};
  phmap::flat_hash_map<ChromosomePair, float> _matrix_tot_counts{};
  ExpectedValuesAggregator _expected_values{};

 public:
  PixelTank() = default;
  explicit PixelTank(std::shared_ptr<const BinTable> bins);

  [[nodiscard]] bool contains(const Chromosome& chrom1, const Chromosome& chrom2) const noexcept;

  void add_pixel(const ThinPixel<N>& p, bool update_expected_values = true);
  void add_pixel(const Pixel<N>& p, bool update_expected_values = true);
  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void add_pixels(PixelIt first_pixel, PixelIt last_pixel, bool update_expected_values = true);

  void finalize();

  [[nodiscard]] auto pixels() const noexcept
      -> const phmap::flat_hash_map<ChromosomePair, ChromPixelTank>&;
  [[nodiscard]] auto matrix_counts() const noexcept
      -> const phmap::flat_hash_map<ChromosomePair, float>&;
  [[nodiscard]] auto expected_values() const noexcept;

  [[nodiscard]] auto pixels(const Chromosome& chrom1, const Chromosome& chrom2) const
      -> const ChromPixelTank&;
  [[nodiscard]] auto matrix_counts(const Chromosome& chrom1, const Chromosome& chrom2) const
      -> const float&;
};

class MetadataOffsetTank {
  struct Key {
    Chromosome chrom1{};
    Chromosome chrom2{};
    std::uint32_t resolution{};

    bool operator<(const Key& other) const noexcept;
  };

  struct Value {
    std::streamoff matrix_metadata_offset{};
    std::size_t matrix_size{};
  };

  phmap::btree_map<Key, Value> _tank{};

 public:
  MetadataOffsetTank() = default;

  [[nodiscard]] bool contains(const Chromosome& chrom1, const Chromosome& chrom2,
                              std::uint32_t resolution) const noexcept;
  [[nodiscard]] auto at(const Chromosome& chrom1, const Chromosome& chrom2,
                        std::uint32_t resolution) const -> const Value&;

  void insert(const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t resolution,
              std::size_t offset, std::size_t size);

  auto operator()() const noexcept -> const phmap::btree_map<Key, Value>&;
};

class HiCBlockPartitioner {
 public:
  class BlockMapperIntra;
  class BlockMapperInter;

 private:
  struct BlockID {
    std::uint32_t chrom1_id;
    std::uint32_t chrom2_id;
    std::uint64_t bid;

    [[nodiscard]] bool operator<(const BlockID& other) const noexcept;
  };

  struct BlockIndex {
    std::uint64_t offset;
    std::uint32_t size;
  };

  std::filesystem::path _path;
  filestream::FileStream _fs;
  std::shared_ptr<const BinTable> _bin_table{};

  phmap::btree_map<BlockID, std::vector<BlockIndex>> _block_index{};
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, std::vector<BlockID>> _chromosome_index{};

  phmap::btree_map<BlockID, MatrixInteractionBlockFlat<float>> _blocks{};
  std::size_t _pixels_processed{};

  phmap::flat_hash_map<Chromosome, BlockMapperIntra> _mappers_intra{};
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, BlockMapperInter> _mappers_inter{};

  BinaryBuffer _bbuffer{};
  int _compression_lvl{};
  std::unique_ptr<ZSTD_CCtx_s> _zstd_cctx{};
  std::unique_ptr<ZSTD_DCtx_s> _zstd_dctx{};
  std::string _compression_buffer{};

  static constexpr std::int32_t DEFAULT_INTRA_CUTOFF = 500;
  static constexpr std::int32_t DEFAULT_INTER_CUTOFF = 5'000;
  static constexpr std::size_t DEFAULT_BLOCK_CAPACITY = 1'000;

 public:
  HiCBlockPartitioner(std::filesystem::path path, std::shared_ptr<const BinTable> bins,
                      int compression_lvl);

  const Reference& chromosomes() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(PixelIt first_pixel, PixelIt last_pixel, std::size_t chunk_size = 100'000'000);

  [[nodiscard]] auto block_index() const noexcept
      -> phmap::btree_map<BlockID, std::vector<BlockIndex>>;
  [[nodiscard]] auto merge_blocks(const BlockID& bid) -> MatrixInteractionBlock<float>;

  void finalize();

 private:
  void init_block_mappers();

  template <typename N>
  [[nodiscard]] auto map(const ThinPixel<N>& p) const -> BlockID;
  template <typename N>
  [[nodiscard]] auto map(const Pixel<N>& p) const -> BlockID;

  void write_blocks();
  void index_chromosomes();
  std::pair<std::uint64_t, std::uint32_t> write_block(const MatrixInteractionBlockFlat<float>& blk);

  [[nodiscard]] std::size_t compute_block_column_count(
      std::size_t num_bins, std::uint32_t bin_size, std::uint32_t cutoff,
      std::size_t block_capacity = DEFAULT_BLOCK_CAPACITY);
  [[nodiscard]] std::size_t compute_num_bins(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                             std::size_t bin_size);

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

struct HiCSectionOffsets {
  std::streamoff position;
  std::size_t size;
};

class ChromChromHiCFileWriter {
  Chromosome _chrom1{};
  Chromosome _chrom2{};

  std::shared_ptr<const HiCHeader> _header{};
  std::shared_ptr<filestream::FileStream> _fs{};
  std::shared_ptr<const BinTable> _bin_table{};
  phmap::btree_map<BlockIndexKey, phmap::btree_set<MatrixBlockMetadata>> _block_index{};
  MatrixMetadata _matrix_metadata{};
  MatrixResolutionMetadata _matrix_resolution_metadata{};
  FooterV5 _footer{};

  BinaryBuffer _bbuffer{};
  std::unique_ptr<libdeflate_compressor> _compressor{};
  std::string _compression_buffer{};

  PixelTank<float> _pixel_tank{};

  std::unique_ptr<const hictk::internal::TmpDir> _tmpdir{};

  HiCSectionOffsets _header_offset{};
  HiCSectionOffsets _body_offset{};
  HiCSectionOffsets _body_metadata_offset{};

 public:
  ChromChromHiCFileWriter() = default;
  explicit ChromChromHiCFileWriter(
      const Chromosome& chrom1, const Chromosome& chrom2, HiCHeader header,
      const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
      std::int32_t compression_lvl = 9, std::size_t buffer_size = 32'000'000);

  [[nodiscard]] std::string_view url() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins() const;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void write_pixels(PixelIt first_pixel, PixelIt last_pixel, bool update_expected_values = true);

  // Write header
  void write_header();
  void write_footer_offset(std::streamoff master_index_offset);
  void write_norm_vector_index(std::streamoff position, std::size_t length);

  // Write body
  auto write_body_metadata(const std::string& unit = "BP") -> HiCSectionOffsets;

  auto write_interaction_block(std::uint64_t block_id, const std::vector<ThinPixel<float>>& pixels,
                               std::size_t bin_column_offset, std::size_t bin_row_offset)
      -> HiCSectionOffsets;

  // Write footer
  auto write_footer() -> HiCSectionOffsets;
  void write_footer_section_size(std::streamoff footer_offset, std::uint64_t bytes);

  // Write expected/normalization values
  void write_expected_values(std::string_view unit);

  // copy sections
  void copy_body(std::ofstream& dest, std::size_t chunk_size = 32'000'000);

  void finalize();

 private:
  [[nodiscard]] static std::shared_ptr<const HiCHeader> init_header(const Chromosome& chrom1,
                                                                    const Chromosome& chrom2,
                                                                    HiCHeader&& header);

  auto write_matrix_metadata() -> HiCSectionOffsets;
  auto write_resolution_metadata(float sum_counts, const std::string& unit) -> HiCSectionOffsets;
};

class HiCFileWriter {
  std::shared_ptr<const HiCHeader> _header{};
  std::shared_ptr<filestream::FileStream> _fs{};
  phmap::flat_hash_map<std::uint32_t, std::shared_ptr<const BinTable>> _bin_tables{};
  phmap::btree_map<BlockIndexKey, phmap::btree_set<MatrixBlockMetadata>> _block_index{};
  std::vector<MatrixMetadata> _matrix_metadata{};
  std::vector<MatrixResolutionMetadata> _matrix_resolution_metadata{};
  std::vector<FooterV5> _footers{};

  BinaryBuffer _bbuffer{};
  std::unique_ptr<libdeflate_compressor> _compressor{};
  std::string _compression_buffer{};

  MetadataOffsetTank _metadata_offset_tank{};
  phmap::flat_hash_map<std::uint32_t, PixelTank<float>> _pixel_tank{};

  std::unique_ptr<const hictk::internal::TmpDir> _tmpdir{};

  struct HiCSectionOffsets {
    std::streamoff position;
    std::size_t size;
  };

  HiCSectionOffsets _header_offsets{};
  HiCSectionOffsets _body_offsets{};
  [[maybe_unused]] HiCSectionOffsets _footer_offsets{};  // TODO

  static constexpr std::int32_t DEFAULT_INTRA_CUTOFF = 500;
  static constexpr std::int32_t DEFAULT_INTER_CUTOFF = 5'000;
  static constexpr std::size_t DEFAULT_BLOCK_CAPACITY = 1'000;

 public:
  HiCFileWriter() = default;
  explicit HiCFileWriter(
      HiCHeader header,
      const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
      std::int32_t compression_lvl = 9, std::size_t buffer_size = 32'000'000);

  [[nodiscard]] std::string_view url() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins(std::uint32_t resolution) const;
  [[nodiscard]] const std::vector<std::uint32_t>& resolutions() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(std::uint32_t resolution, PixelIt first_pixel, PixelIt last_pixel,
                     bool update_expected_values = true);

  void write_pixels(bool write_chromosome_ALL = true);
  void write_pixels(const Chromosome& chrom1, const Chromosome& chrom2);
  void write_pixels(const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t resolution);
  void write_pixels_ALL(std::size_t num_bins = 1000);

  // Write header
  void write_header();
  void write_footer_offset(std::streamoff master_index_offset);
  void write_norm_vector_index(std::streamoff position, std::size_t length);

  // Write body
  auto write_body_metadata(std::uint32_t resolution, const Chromosome& chrom1,
                           const Chromosome& chrom2, const std::string& unit = "BP")
      -> HiCSectionOffsets;

  auto write_interaction_block(std::uint64_t block_id, const Chromosome& chrom1,
                               const Chromosome& chrom2, std::uint32_t resolution,
                               const std::vector<ThinPixel<float>>& pixels,
                               std::size_t bin_column_offset, std::size_t bin_row_offset)
      -> HiCSectionOffsets;

  // Write footer
  auto write_footers() -> HiCSectionOffsets;
  void add_footer(const Chromosome& chrom1, const Chromosome& chrom2, std::size_t file_offset,
                  std::size_t matrix_metadata_bytes);
  void write_footer_section_size(std::streamoff footer_offset, std::uint64_t bytes);

  // Write expected/normalization values
  void write_expected_values(std::string_view unit);

  // copy sections
  void copy_body(HiCFileWriter& dest, std::size_t chunk_size = 32'000'000);

  void finalize();

 private:
  [[nodiscard]] std::size_t compute_block_column_count(
      std::size_t num_bins, std::uint32_t bin_size, std::uint32_t cutoff,
      std::size_t block_capacity = DEFAULT_BLOCK_CAPACITY);
  [[nodiscard]] std::size_t compute_num_bins(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                             std::size_t bin_size);

  auto write_matrix_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id) -> HiCSectionOffsets;
  auto write_resolutions_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                  float sum_counts, const std::string& unit) -> HiCSectionOffsets;

  void coarsen_pixels();
};
}  // namespace hictk::hic::internal

#include "./impl/file_writer_impl.hpp"  // NOLINT
