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
  using Row = std::vector<Pixel<N>>;

 public:
  std::int32_t nRecords{};
  std::int32_t binColumnOffset{std::numeric_limits<std::int32_t>::max()};
  std::int32_t binRowOffset{std::numeric_limits<std::int32_t>::max()};
  std::uint8_t useFloatContact{};
  std::uint8_t useIntXPos{};
  std::uint8_t useIntYPos{};
  std::uint8_t matrixRepresentation{};

  void emplace_back(Pixel<N>&& p);
  void finalize();

  [[nodiscard]] auto operator()() const noexcept -> const phmap::btree_map<RowID, Row>&;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, libdeflate_compressor& compressor,
                                      std::string& compression_buffer, bool clear = true) const;

 private:
  phmap::btree_map<RowID, Row> _interactions;
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

}  // namespace hictk::hic::internal

#include "./impl/file_writer_data_structures_impl.hpp"  // NOLINT
