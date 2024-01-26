// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"
#include <libdeflate.h>
#include <parallel_hashmap/btree.h>

#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/hic/binary_buffer.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/pixel.hpp"

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

  [[nodiscard]] bool operator<(const MatrixResolutionMetadata& other) const noexcept;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;

  template <typename It>
  void set_block_metadata(It first_block, It last_block);

 private:
  std::vector<MatrixBlockMetadata> _block_metadata{};
};

struct MatrixBodyMetadata {
  MatrixMetadata matrixMetadata;
  phmap::btree_set<MatrixResolutionMetadata> resolutionMetadata;

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#blocks
template <typename N = float>
struct MatrixInteractionBlock {
 private:
  using RowID = std::int32_t;

  struct Pixel {
    std::int32_t column;
    N count;
    [[nodiscard]] bool operator<(const Pixel& other) const noexcept;
  };
  using Row = phmap::btree_set<Pixel>;

 public:
  std::int32_t nRecords{};
  std::int32_t binColumnOffset{std::numeric_limits<std::int32_t>::max()};
  std::int32_t binRowOffset{std::numeric_limits<std::int32_t>::max()};
  std::uint8_t useFloatContact{};
  std::uint8_t useIntXPos{};
  std::uint8_t useIntYPos{};
  std::uint8_t matrixRepresentation{};

  std::int16_t w{};

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] double sum() const noexcept;

  void emplace_back(hictk::Pixel<N>&& p, std::uint32_t bin_id_offset = 0);
  void finalize();

  [[nodiscard]] auto operator()() const noexcept -> const phmap::btree_map<RowID, Row>&;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, libdeflate_compressor& compressor,
                                      std::string& compression_buffer, bool clear = true) const;

 private:
  double _sum{};
  phmap::btree_map<RowID, Row> _interactions;

  std::int32_t _min_col{std::numeric_limits<std::int32_t>::max()};
  std::int32_t _max_col{};

  [[nodiscard]] std::size_t compute_size_lor_repr() const noexcept;
  [[nodiscard]] std::size_t compute_size_dense_repr() const noexcept;

  [[nodiscard]] std::size_t compute_dense_width() const noexcept;

  [[nodiscard]] std::string serialize_lor(BinaryBuffer& buffer, libdeflate_compressor& compressor,
                                          std::string& compression_buffer, bool clear = true) const;
  [[nodiscard]] std::string serialize_dense(BinaryBuffer& buffer, libdeflate_compressor& compressor,
                                            std::string& compression_buffer,
                                            bool clear = true) const;

  static void compress(const std::string& buffer_in, std::string& buffer_out,
                       libdeflate_compressor& compressor);
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#master-index
struct FooterMasterIndex {
  std::string key;
  std::int64_t position;
  std::int32_t size;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

struct ExpectedValuesBlock {
  std::string unit{};
  std::int32_t binSize{};
  [[nodiscard]] std::int64_t nValues() const noexcept;
  std::vector<float> value{};
  [[nodiscard]] std::int32_t nChrScaleFactors() const noexcept;
  std::vector<std::int32_t> chrIndex{};
  std::vector<float> chrScaleFactor{};

  ExpectedValuesBlock() = default;
  ExpectedValuesBlock(std::string_view unit_, std::uint32_t bin_size,
                      const std::vector<double>& weights,
                      const std::vector<std::uint32_t>& chrom_ids,
                      const std::vector<double>& scale_factors);

  [[nodiscard]] bool operator<(const ExpectedValuesBlock& other) const noexcept;

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
  [[nodiscard]] static ExpectedValuesBlock deserialize(filestream::FileStream& fs);
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#expected-value-vectors
class ExpectedValues {
  phmap::btree_set<ExpectedValuesBlock> _expected_values;

 public:
  [[nodiscard]] std::int32_t nExpectedValueVectors() const noexcept;
  [[nodiscard]] const phmap::btree_set<ExpectedValuesBlock>& expectedValues() const noexcept;
  void emplace(const ExpectedValuesBlock& evb, bool force_overwrite = false);
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
  [[nodiscard]] static ExpectedValues deserialize(filestream::FileStream& fs);
};

struct NormalizedExpectedValuesBlock {
  std::string type{};
  std::string unit{};
  std::int32_t binSize{};
  [[nodiscard]] std::int64_t nValues() const noexcept;
  std::vector<float> value{};
  [[nodiscard]] std::int32_t nChrScaleFactors() const noexcept;
  std::vector<std::int32_t> chrIndex{};
  std::vector<float> chrScaleFactor{};

  NormalizedExpectedValuesBlock() = default;
  NormalizedExpectedValuesBlock(std::string_view type_, std::string_view unit_,
                                std::uint32_t bin_size, const std::vector<double>& weights,
                                const std::vector<std::uint32_t>& chrom_ids,
                                const std::vector<double>& scale_factors);

  [[nodiscard]] bool operator<(const NormalizedExpectedValuesBlock& other) const noexcept;

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
  [[nodiscard]] static NormalizedExpectedValuesBlock deserialize(filestream::FileStream& fs);
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalized-expected-value-vectors
class NormalizedExpectedValues {
  phmap::btree_set<NormalizedExpectedValuesBlock> _normalized_expected_values;

 public:
  [[nodiscard]] std::int32_t nNormExpectedValueVectors() const noexcept;
  [[nodiscard]] const phmap::btree_set<NormalizedExpectedValuesBlock>& normExpectedValues()
      const noexcept;
  void emplace(const NormalizedExpectedValuesBlock& evb, bool force_overwrite = false);
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
  [[nodiscard]] static NormalizedExpectedValues deserialize(filestream::FileStream& fs);
};

struct NormalizationVectorIndexBlock {
  std::string type{};
  std::int32_t chrIdx{};
  std::string unit{};
  std::int32_t binSize{};
  std::int64_t position{};
  std::int64_t nBytes{};

 private:
 public:
  NormalizationVectorIndexBlock() = default;
  NormalizationVectorIndexBlock(std::string type_, std::uint32_t chrom_idx, std::string unit_,
                                std::uint32_t bin_size, std::size_t position_, std::size_t n_bytes);

  [[nodiscard]] bool operator<(const NormalizationVectorIndexBlock& other) const noexcept;

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
  [[nodiscard]] static NormalizationVectorIndexBlock deserialize(filestream::FileStream& fs);
};

// https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalization-vector-index
class NormalizationVectorIndex {
  std::vector<NormalizationVectorIndexBlock> _norm_vect_idx{};

 public:
  [[nodiscard]] std::int32_t nNormVectors() const noexcept;
  [[nodiscard]] const std::vector<NormalizationVectorIndexBlock> normalizationVectorIndex()
      const noexcept;
  void emplace_back(NormalizationVectorIndexBlock blk);

  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
  [[nodiscard]] static NormalizationVectorIndex deserialize(filestream::FileStream& fs);
};

}  // namespace hictk::hic::internal

#include "./impl/file_writer_data_structures_impl.hpp"  // NOLINT
