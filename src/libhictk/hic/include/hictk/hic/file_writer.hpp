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

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
};

struct MatrixBlockMetadata {
  std::int32_t blockNumber{};
  std::int64_t blockPosition{};
  std::int32_t blockSizeBytes{};

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
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

  std::vector<MatrixBlockMetadata> blocksMetadata{};

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
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

  MatrixInteractionBlock(const InteractionBlock& blk, std::size_t bin_row_offset);

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, libdeflate_compressor& compressor,
                                      std::string& compression_buffer) const;

 private:
  using RowID = std::int32_t;
  using Row = std::vector<ThinPixel<float>>;
  phmap::btree_map<RowID, Row> _interactions;

  auto group_interactions_by_column(const InteractionBlock& blk, std::size_t bin_row_offset)
      -> phmap::btree_map<RowID, Row>;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#master-index
struct MasterIndex {
  std::string key;
  std::int64_t position;
  std::int32_t size;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#expected-value-vectors
struct ExpectedValues {
  std::int32_t nExpectedValueVectors = 0;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalized-expected-value-vectors
struct NormalizedExpectedValues {
  std::int32_t nNormExpectedValueVectors = 0;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalization-vector-index
struct NormalizationVectorIndex {
  std::int32_t nNormVectors = 0;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalization-vector-arrays-1-per-normalization-vector
struct NormalizationVectorArray {
  std::int64_t nValues = 0;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
};

struct FooterV5 {
  std::int64_t nBytesV5;
  std::int32_t nEntries;

  std::vector<MasterIndex> masterIndex;

  ExpectedValues expectedValues;
  NormalizedExpectedValues normExpectedValues;
  NormalizationVectorIndex normVectIndex;
  std::vector<NormalizationVectorArray> normVectArray;

  FooterV5() = default;
  [[nodiscard]] constexpr std::int64_t master_index_offset() const noexcept;
  void add_footer(const HiCFooter& footer, std::int64_t matrix_offset, std::int32_t matrix_bytes);
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer) const;
};

class HiCFileWriter {
  std::shared_ptr<const HiCHeader> _header{};
  std::shared_ptr<filestream::FileStream> _fs{};
  phmap::btree_set<MatrixBlockMetadata> _block_index{};

  BinaryBuffer _bbuffer{};
  std::unique_ptr<libdeflate_compressor> _compressor{};
  std::string _compression_buffer{};

  static constexpr std::int32_t DEFAULT_INTRA_CUTOFF = 500;
  static constexpr std::int32_t DEFAULT_INTER_CUTOFF = 5'000;
  static constexpr std::size_t DEFAULT_BLOCK_CAPACITY = 1'000;

 public:
  HiCFileWriter() = default;
  explicit HiCFileWriter(HiCHeader header, std::int32_t compression_lvl = 9,
                         std::size_t buffer_size = 32'000'000);

  [[nodiscard]] std::string_view url() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const std::vector<std::uint32_t> resolutions() const noexcept;

  // Write header
  void write_header();
  void write_master_index_offset(std::int64_t master_index);

  // Write body
  std::pair<std::streamoff, std::size_t> write_body_metadata(std::uint32_t chrom1_id,
                                                             std::uint32_t chrom2_id,
                                                             const std::string& unit = "BP");

  std::streamoff write_interaction_block(const InteractionBlock& blk, std::size_t bin_column_offset,
                                         std::size_t bin_row_offset);

  std::streamoff write_footer(const std::vector<HiCFooter>& footers,
                              const std::vector<std::int64_t>& matrix_offsets,
                              const std::vector<std::int32_t>& matrix_bytes);

 private:
  [[nodiscard]] std::size_t compute_block_column_count(
      std::size_t num_bins, std::uint32_t bin_size, std::uint32_t cutoff,
      std::size_t block_capacity = DEFAULT_BLOCK_CAPACITY);
  [[nodiscard]] std::size_t compute_num_bins(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                             std::size_t bin_size);
  [[nodiscard]] static phmap::btree_map<std::int32_t, std::vector<ThinPixel<float>>>
  group_interactions_by_column(const InteractionBlock& blk, std::size_t bin_row_offset);

  void write_matrix_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id);
  void write_resolution_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                 const std::string& unit);

 public:
  class BlockMapperInter {
    std::uint64_t _block_bin_count{};
    std::uint64_t _block_column_count{};

   public:
    BlockMapperInter(std::uint64_t block_bin_count, std::uint64_t block_column_count);
    [[nodiscard]] std::uint64_t operator()(std::uint64_t bin1_id, std::uint64_t bin2_id);

    [[nodiscard]] std::uint64_t block_bin_count();
    [[nodiscard]] std::uint64_t block_column_count();
  };

  class BlockMapperIntra {
    BlockMapperInter _inter_mapper;
    double _base{};

    static constexpr std::int64_t DEFAULT_BASE_DEPTH = 2;

   public:
    BlockMapperIntra(std::uint64_t block_bin_count, std::uint64_t block_column_count,
                     std::int64_t base_depth = DEFAULT_BASE_DEPTH);
    [[nodiscard]] std::uint64_t operator()(std::uint64_t bin1_id, std::uint64_t bin2_id);

    [[nodiscard]] std::uint64_t block_bin_count();
    [[nodiscard]] std::uint64_t block_column_count();

   private:
    [[nodiscard]] bool use_inter_mapper() const noexcept;
    [[nodiscard]] static double init_base(std::int64_t base_depth) noexcept;
  };
};
}  // namespace hictk::hic::internal

#include "./impl/file_writer_impl.hpp"  // NOLINT
