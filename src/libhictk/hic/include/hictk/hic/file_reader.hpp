// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <libdeflate.h>

#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/hic/header.hpp"
#include "hictk/hic/index.hpp"

namespace hictk::hic::internal {

class HiCFileReader {
  using Decompressor = UniquePtrWithDeleter<libdeflate_decompressor>;
  std::shared_ptr<filestream::FileStream> _fs{};
  std::shared_ptr<const HiCHeader> _header{};
  std::string _strbuff{};
  Decompressor _decompressor{init_decompressor()};

 public:
  HiCFileReader() = default;
  explicit HiCFileReader(std::string url);
  [[nodiscard]] inline const std::string &url() const noexcept;
  [[nodiscard]] const HiCHeader &header() const noexcept;

  [[nodiscard]] std::int32_t version() const noexcept;

  // reads the footer given a pair of chromosomes, wanted_norm, wanted_unit (BP or FRAG) and
  // resolution.
  [[nodiscard]] HiCFooter read_footer(
      std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
      balancing::Method wanted_norm, MatrixUnit wanted_unit, std::uint32_t wanted_resolution,
      std::shared_ptr<balancing::Weights> weights1 = std::make_shared<balancing::Weights>(),
      std::shared_ptr<balancing::Weights> weights2 = std::make_shared<balancing::Weights>());

  [[nodiscard]] std::int64_t read_footer_file_offset(std::string_view key);
  [[nodiscard]] std::vector<double> read_footer_expected_values(
      std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
      balancing::Method wanted_norm, MatrixUnit wanted_unit, std::uint32_t wanted_resolution);
  [[nodiscard]] std::vector<double> read_footer_expected_values_norm(
      std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
      balancing::Method wanted_norm, MatrixUnit wanted_unit, std::uint32_t wanted_resolution);
  void read_footer_norm(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                        balancing::Method wanted_norm, MatrixUnit wanted_unit,
                        std::uint32_t wanted_resolution, const Chromosome &chrom1,
                        const Chromosome &chrom2, std::shared_ptr<balancing::Weights> &weights1,
                        std::shared_ptr<balancing::Weights> &weights2);

  [[nodiscard]] std::vector<balancing::Method> list_avail_normalizations(
      MatrixType matrix_type, MatrixUnit wanted_unit, std::uint32_t wanted_resolution);

  [[nodiscard]] static MatrixType readMatrixType(filestream::FileStream &fs, std::string &buff);
  [[nodiscard]] static balancing::Method readNormalizationMethod(filestream::FileStream &fs,
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
  [[nodiscard]] balancing::Method readNormalizationMethod();
  [[nodiscard]] MatrixUnit readMatrixUnit();

  [[nodiscard]] std::int64_t readNValues();
  [[nodiscard]] bool checkMagicString();
  [[nodiscard]] static bool checkMagicString(filestream::FileStream &fs);
  [[nodiscard]] std::int64_t masterOffset() const noexcept;

  [[nodiscard]] static auto init_decompressor() -> Decompressor;
};
}  // namespace hictk::hic::internal

#include "./impl/file_reader_impl.hpp"
