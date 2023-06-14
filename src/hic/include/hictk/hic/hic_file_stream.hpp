// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <libdeflate.h>
#include <parallel_hashmap/btree.h>

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/hic_footer.hpp"
#include "hictk/hic/hic_header.hpp"
#include "hictk/hic/index.hpp"

namespace hictk::hic::internal {

class HiCFileStream {
  using Decompressor = UniquePtrWithDeleter<libdeflate_decompressor>;
  std::shared_ptr<filestream::FileStream> _fs{};
  std::shared_ptr<const HiCHeader> _header{};
  std::string _strbuff{};
  Decompressor _decompressor{init_decompressor()};

 public:
  HiCFileStream() = default;
  explicit HiCFileStream(std::string url);
  [[nodiscard]] inline const std::string &url() const noexcept;
  [[nodiscard]] const HiCHeader &header() const noexcept;

  [[nodiscard]] std::int32_t version() const noexcept;

  // reads the footer given a pair of chromosomes, wanted_norm, wanted_unit (BP or FRAG) and
  // resolution.
  [[nodiscard]] HiCFooter read_footer(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                      MatrixType matrix_type, NormalizationMethod wanted_norm,
                                      MatrixUnit wanted_unit, std::uint32_t wanted_resolution);

  [[nodiscard]] static MatrixType readMatrixType(filestream::FileStream &fs, std::string &buff);
  [[nodiscard]] static NormalizationMethod readNormalizationMethod(filestream::FileStream &fs,
                                                                   std::string &buff);
  [[nodiscard]] static MatrixUnit readMatrixUnit(filestream::FileStream &fs, std::string &buff);

  [[nodiscard]] Index read_index(std::int64_t fileOffset, const Chromosome &chrom1,
                                 const Chromosome &chrom2, MatrixUnit wantedUnit,
                                 std::int64_t wantedResolution);
  void readAndInflate(const BlockIndex &idx, std::string &plainTextBuffer);

  [[nodiscard]] static bool checkMagicString(std::string url) noexcept;

 private:
  [[nodiscard]] static filestream::FileStream openStream(std::string url);
  // reads the header, storing the positions of the normalization vectors and returning the
  // masterIndexPosition pointer
  [[nodiscard]] static HiCHeader readHeader(filestream::FileStream &fs);

  [[nodiscard]] std::vector<double> readExpectedVector(std::int64_t nValues);
  [[nodiscard]] std::vector<double> readNormalizationFactors(std::uint32_t wantedChrom);
  void applyNormalizationFactors(std::vector<double> &expectedValues,
                                 const std::vector<double> &normFactors);
  [[nodiscard]] std::vector<double> readNormalizationVector(indexEntry cNormEntry,
                                                            std::size_t numValuesExpected);

  void discardExpectedVector(std::int64_t nValues);
  void discardNormalizationFactors(std::uint32_t wantedChrom);

  [[nodiscard]] MatrixType readMatrixType();
  [[nodiscard]] NormalizationMethod readNormalizationMethod();
  [[nodiscard]] MatrixUnit readMatrixUnit();

  [[nodiscard]] std::int64_t readNValues();
  [[nodiscard]] bool checkMagicString();
  [[nodiscard]] static bool checkMagicString(filestream::FileStream &fs);
  [[nodiscard]] std::int64_t masterOffset() const noexcept;

  [[nodiscard]] static auto init_decompressor() -> Decompressor;
};
}  // namespace hictk::hic::internal

#include "../../../hic_file_stream_impl.hpp"
