// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <libdeflate.h>
#include <parallel_hashmap/btree.h>

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

template <>
struct std::default_delete<libdeflate_compressor> {
  void operator()(libdeflate_compressor* compressor) const {
    libdeflate_free_compressor(compressor);
  }
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
struct MatrixInteractionBlock {
  std::int32_t nRecords;
  std::int32_t binColumnOffset;
  std::int32_t binRowOffset;
  std::uint8_t useFloatContact;
  std::uint8_t useIntXPos;
  std::uint8_t useIntYPos;
  std::uint8_t matrixRepresentation;

  MatrixInteractionBlock(const BinTable& bins, const std::vector<ThinPixel<float>>& pixels,
                         std::size_t bin_row_offset);

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, libdeflate_compressor& compressor,
                                      std::string& compression_buffer, bool clear = true) const;

 private:
  using RowID = std::int32_t;
  using Row = std::vector<Pixel<float>>;
  phmap::btree_map<RowID, Row> _interactions;

  auto group_interactions_by_column(const BinTable& bins,
                                    const std::vector<ThinPixel<float>>& pixels)
      -> phmap::btree_map<RowID, Row>;
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

  using PixelTankKey = std::pair<Chromosome, Chromosome>;
  using ChromPixelTank = phmap::btree_set<ThinPixel<float>>;
  using PixelTank = phmap::flat_hash_map<PixelTankKey, ChromPixelTank>;
  PixelTank _pixel_tank{};

  phmap::flat_hash_map<std::uint32_t, ExpectedValuesAggregator> _expected_values{};

  static constexpr std::int32_t DEFAULT_INTRA_CUTOFF = 500;
  static constexpr std::int32_t DEFAULT_INTER_CUTOFF = 5'000;
  static constexpr std::size_t DEFAULT_BLOCK_CAPACITY = 1'000;

 public:
  HiCFileWriter() = default;
  explicit HiCFileWriter(HiCHeader header, std::int32_t compression_lvl = 9,
                         std::size_t buffer_size = 32'000'000);

  [[nodiscard]] std::string_view url() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins(std::uint32_t resolution) const;
  [[nodiscard]] const std::vector<std::uint32_t> resolutions() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(std::uint32_t resolution, PixelIt first_pixel, PixelIt last_pixel,
                     bool update_expected_values = true);

  void write_pixels(bool write_chromosome_ALL = true);
  void write_pixels(const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t resolution);
  void write_pixels_ALL(std::size_t num_bins = 1000);

  // Write header
  void write_header();
  void write_footer_offset(std::int64_t master_index_offset);
  void write_norm_vector_index(std::size_t position, std::size_t length);

  // Write body
  auto write_body_metadata(const Chromosome& chrom1, const Chromosome& chrom2,
                           const std::string& unit = "BP");

  std::streamoff write_interaction_block(std::uint64_t block_id, const Chromosome& chrom1,
                                         const Chromosome& chrom2, std::uint32_t resolution,
                                         const std::vector<ThinPixel<float>>& pixels,
                                         std::size_t bin_column_offset, std::size_t bin_row_offset);

  // Write footer
  std::streamoff write_footers();
  void add_footer(const Chromosome& chrom1, const Chromosome& chrom2, std::size_t file_offset,
                  std::size_t matrix_metadata_bytes);
  void write_footer_section_size(std::uint64_t footer_offset, std::uint64_t bytes);

  // Write expected/normalization values
  void write_expected_values(std::string_view unit);

  void finalize();

 private:
  [[nodiscard]] std::size_t compute_block_column_count(
      std::size_t num_bins, std::uint32_t bin_size, std::uint32_t cutoff,
      std::size_t block_capacity = DEFAULT_BLOCK_CAPACITY);
  [[nodiscard]] std::size_t compute_num_bins(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                             std::size_t bin_size);

  std::size_t write_matrix_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id);
  auto write_resolutions_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                  const std::string& unit);

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

#include "./impl/file_writer_impl.hpp"  // NOLINT
