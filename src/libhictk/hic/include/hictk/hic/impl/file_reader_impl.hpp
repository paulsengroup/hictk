// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/compile.h>
#include <fmt/format.h>
#include <libdeflate.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/filestream.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/reference.hpp"

namespace hictk::hic::internal {

inline HiCFileReader::HiCFileReader(std::string url)
    : _fs(std::make_shared<filestream::FileStream<>>(HiCFileReader::openStream(std::move(url)))),
      _header(std::make_shared<const HiCHeader>(HiCFileReader::readHeader(*_fs))) {}

inline filestream::FileStream<> HiCFileReader::openStream(std::string url) {
  try {
    return {std::move(url), nullptr};  // opening wo/ locking is ok
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Failed to open file: {}"), e.what()));
  }
}

inline const std::string &HiCFileReader::path() const noexcept { return _fs->path(); }
inline const HiCHeader &HiCFileReader::header() const noexcept { return *_header; }

inline std::int32_t HiCFileReader::version() const noexcept {
  assert(_header->version != -1);
  return _header->version;
}

inline void HiCFileReader::discardExpectedVector(std::int64_t nValues) {
  const std::int64_t elementSize = version() > 8 ? sizeof(float) : sizeof(double);
  _fs->seekg(nValues * elementSize, std::ios::cur);
}

inline std::vector<double> HiCFileReader::readExpectedVector(std::int64_t nValues) {
  std::vector<double> initialExpectedValues(static_cast<std::size_t>(nValues));
  if (version() > 8) {  // NOLINT(*-avoid-magic-numbers)
    std::vector<float> tmpbuff(static_cast<std::size_t>(nValues));
    _fs->read(tmpbuff);
    std::transform(tmpbuff.begin(), tmpbuff.end(), initialExpectedValues.begin(),
                   [](float n) { return static_cast<double>(n); });
  } else {
    _fs->read(initialExpectedValues);
  }

  // This seems to be copying initialValues into finalResult at the moment
  // std::int32_t window = 5000000 / resolution;
  // rollingMedian(initialExpectedValues, _expectedValues, window);
  return initialExpectedValues;
}

inline std::vector<double> HiCFileReader::readNormalizationFactors(std::uint32_t wantedChrom) {
  const auto nFactors = _fs->read<std::int32_t>();
  std::vector<double> normFactors{};
  auto readFactor = [this]() {
    if (version() > 8) {  // NOLINT(*-avoid-magic-numbers)
      return _fs->read_as_double<float>();
    }
    return _fs->read<double>();
  };

  for (auto i = 0; i < nFactors; ++i) {
    const auto foundChrom = _fs->read_as_unsigned<std::int32_t>();
    const auto v = readFactor();
    if (foundChrom == wantedChrom) {
      normFactors.push_back(v);
    }
  }
  return normFactors;
}

inline void HiCFileReader::applyNormalizationFactors(std::vector<double> &expectedValues,
                                                     const std::vector<double> &normFactors) {
  if (normFactors.empty() || expectedValues.empty()) {
    return;
  }
  for (const auto factor : normFactors) {
    std::transform(expectedValues.begin(), expectedValues.end(), expectedValues.begin(),
                   [&](auto n) { return n / factor; });
  }
}
inline std::vector<double> HiCFileReader::readNormalizationVector(indexEntry cNormEntry,
                                                                  std::size_t numValuesExpected) {
  _fs->seekg(cNormEntry.position);
  const auto numValues = static_cast<std::size_t>(readNValues());

  // We cannot use numValues directly because sometimes hic files have few trailing zeros for some
  // reason
  if (numValues < numValuesExpected) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("normalization vector is corrupted: expected {} values, found {}"),
                    numValuesExpected, numValues));
  }

  std::vector<double> buffer(numValuesExpected);
  if (version() > 8) {  // NOLINT(*-avoid-magic-numbers)
    std::vector<float> tmpbuffer(numValuesExpected);
    _fs->read(tmpbuffer);
    std::transform(tmpbuffer.begin(), tmpbuffer.end(), buffer.begin(),
                   [](float n) { return static_cast<double>(n); });

  } else {
    _fs->read(buffer);
  }

  return buffer;
}

inline void HiCFileReader::discardNormalizationFactors(std::uint32_t wantedChrom) {
  std::ignore = readNormalizationFactors(wantedChrom);
}

inline MatrixType HiCFileReader::readMatrixType(filestream::FileStream<> &fs, std::string &buff) {
  fs.getline(buff, '\0');
  return ParseMatrixTypeStr(buff);
}

inline balancing::Method HiCFileReader::readNormalizationMethod(filestream::FileStream<> &fs,
                                                                std::string &buff) {
  fs.getline(buff, '\0');
  return balancing::Method{buff};
}

inline MatrixUnit HiCFileReader::readMatrixUnit(filestream::FileStream<> &fs, std::string &buff) {
  fs.getline(buff, '\0');
  return ParseUnitStr(buff);
}

inline MatrixType HiCFileReader::readMatrixType() { return readMatrixType(*_fs, _strbuff); }

inline balancing::Method HiCFileReader::readNormalizationMethod() {
  return readNormalizationMethod(*_fs, _strbuff);
}

inline MatrixUnit HiCFileReader::readMatrixUnit() { return readMatrixUnit(*_fs, _strbuff); }

inline std::int64_t HiCFileReader::readNValues() {
  if (version() > 8) {  // NOLINT(*-avoid-magic-numbers)
    return _fs->read<std::int64_t>();
  }
  return _fs->read<std::int32_t>();
}

inline bool HiCFileReader::checkMagicString(filestream::FileStream<> &fs) {
  return fs.getline('\0') == "HIC";
}

inline std::int64_t HiCFileReader::masterOffset() const noexcept { return _header->footerPosition; }

inline auto HiCFileReader::init_decompressor() -> Decompressor {
  Decompressor zs(libdeflate_alloc_decompressor(),
                  [](auto *ptr) { libdeflate_free_decompressor(ptr); });
  if (!zs) {
    throw std::runtime_error("failed to initialize zlib decompression stream");
  }

  return zs;
}

inline Index HiCFileReader::read_index(std::int64_t fileOffset, const Chromosome &chrom1,
                                       const Chromosome &chrom2, MatrixUnit wantedUnit,
                                       std::int64_t wantedResolution) {
  _fs->seekg(fileOffset);

  [[maybe_unused]] const auto c1i = _fs->read<std::int32_t>();
  [[maybe_unused]] const auto c2i = _fs->read<std::int32_t>();
  const auto numResolutions = _fs->read<std::int32_t>();

  assert(c1i == static_cast<std::int32_t>(chrom1.id()));
  assert(c2i == static_cast<std::int32_t>(chrom2.id()));

  for (std::int32_t i = 0; i < numResolutions; ++i) {
    const auto foundUnit = readMatrixUnit();
    std::ignore = _fs->read<std::int32_t>();  // oldIndex
    const auto sumCount = _fs->read<float>();
    std::ignore = _fs->read<float>();  // occupiedCellCount
    std::ignore = _fs->read<float>();  // percent5
    std::ignore = _fs->read<float>();  // percent95

    const auto foundResolution = static_cast<std::int64_t>(_fs->read<std::int32_t>());
    const auto blockBinCount = static_cast<std::size_t>(_fs->read<std::int32_t>());
    const auto blockColumnCount = static_cast<std::size_t>(_fs->read<std::int32_t>());

    const auto nBlocks = static_cast<std::size_t>(_fs->read<std::int32_t>());

    if (wantedUnit == foundUnit && wantedResolution == foundResolution) {
      Index::BlkIdxBuffer buffer(nBlocks);
      for (std::size_t j = 0; j < nBlocks; ++j) {
        const auto block_id = static_cast<std::size_t>(_fs->read<std::int32_t>());
        const auto position = static_cast<std::size_t>(_fs->read<std::int64_t>());
        const auto size = static_cast<std::size_t>(_fs->read<std::int32_t>());
        assert(static_cast<std::streamsize>(position + size) < _fs->size());
        if (size > 0) {
          buffer.emplace(block_id, position, size, blockColumnCount);
        }
      }

      return {chrom1,           chrom2,
              wantedUnit,       static_cast<std::uint32_t>(wantedResolution),
              version(),        blockBinCount,
              blockColumnCount, static_cast<double>(sumCount),
              std::move(buffer)};
    }

    constexpr std::int64_t blockSize = sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t);
    _fs->seekg(static_cast<std::int64_t>(nBlocks) * blockSize, std::ios::cur);
  }

  throw std::runtime_error(
      fmt::format(FMT_STRING("Unable to find block map for {}:{} with unit {} and resolution {}"),
                  chrom1.name(), chrom2.name(), wantedUnit, wantedResolution));
}

inline bool HiCFileReader::checkMagicString() { return checkMagicString(*_fs); }

inline HiCHeader HiCFileReader::readHeader(filestream::FileStream<> &fs) {
  return HiCHeader::deserialize(0, fs);
}

inline void HiCFileReader::readAndInflate(const BlockIndex &idx, std::string &plainTextBuffer) {
  try {
    // _strbuff is used to store compressed data
    // plainTextBuffer is used to store decompressed data
    assert(_decompressor);
    assert(idx.compressed_size_bytes() > 0);
    const auto buffSize = idx.compressed_size_bytes();

    plainTextBuffer.reserve(buffSize * 3);
    plainTextBuffer.resize(plainTextBuffer.capacity());

    _fs->seekg(static_cast<std::int64_t>(idx.file_offset()));
    _fs->read(_strbuff, buffSize);

    std::size_t bytes_decompressed{};

    while (true) {
      using LR = libdeflate_result;
      const auto status = libdeflate_zlib_decompress(_decompressor.get(), _strbuff.data(),
                                                     _strbuff.size(), plainTextBuffer.data(),
                                                     plainTextBuffer.size(), &bytes_decompressed);
      if (status == LR::LIBDEFLATE_SUCCESS) {
        plainTextBuffer.resize(bytes_decompressed);
        break;
      }
      if (status == LR::LIBDEFLATE_INSUFFICIENT_SPACE) {
        plainTextBuffer.resize(plainTextBuffer.size() + buffSize);
        continue;
      }
      if (status == LR::LIBDEFLATE_BAD_DATA) {
        throw std::runtime_error("invalid or corrupted data");
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to decompress block at pos {}: {}"),
                                         idx.file_offset(), e.what()));
  }
}

inline bool HiCFileReader::checkMagicString(std::string url) noexcept {
  try {
    filestream::FileStream fs(openStream(std::move(url)));
    return checkMagicString(fs);
  } catch (...) {
    return false;
  }
}

inline std::int64_t HiCFileReader::read_footer_file_offset(std::string_view key) {
  std::ignore = readNValues();  // nBytes

  std::int64_t pos = -1;
  auto nEntries = _fs->read<std::int32_t>();
  for (int i = 0; i < nEntries; i++) {
    const auto strbuff = _fs->getline('\0');
    assert(!strbuff.empty());
    const auto fpos = _fs->read<std::int64_t>();
    std::ignore = _fs->read<std::int32_t>();  // sizeInBytes
    if (strbuff == key) {
      pos = fpos;
    }
  }
  return pos;
}

inline std::vector<double> HiCFileReader::read_footer_expected_values(
    const Chromosome &chrom1, const Chromosome &chrom2, MatrixType matrix_type,
    const balancing::Method &wanted_norm, MatrixUnit wanted_unit, std::uint32_t wanted_resolution) {
  std::vector<double> expectedValues{};
  auto nExpectedValues = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nExpectedValues; ++i) {
    const auto foundUnit = readMatrixUnit();
    const auto foundResolution = _fs->read_as_unsigned<std::int32_t>();
    const auto nValues = readNValues();

    using MT = MatrixType;
    using NM = balancing::Method;
    const bool store = chrom1 == chrom2 && (matrix_type == MT::oe || matrix_type == MT::expected) &&
                       wanted_norm == NM::NONE() && foundUnit == wanted_unit &&
                       foundResolution == wanted_resolution;

    if (store) {
      expectedValues = readExpectedVector(nValues);
      const auto normFactors = readNormalizationFactors(chrom1.id());
      applyNormalizationFactors(expectedValues, normFactors);

    } else {
      discardExpectedVector(nValues);
      discardNormalizationFactors(chrom1.id());
    }
  }
  return expectedValues;
}

inline std::vector<double> HiCFileReader::read_footer_expected_values_norm(
    const Chromosome &chrom1, const Chromosome &chrom2, MatrixType matrix_type,
    const balancing::Method &wanted_norm, MatrixUnit wanted_unit, std::uint32_t wanted_resolution) {
  if (_fs->tellg() == static_cast<std::streampos>(_fs->size())) {
    return {};
  }

  std::vector<double> expectedValues{};
  const auto nExpectedValues = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nExpectedValues; i++) {
    const auto foundNorm = readNormalizationMethod();
    const auto foundUnit = readMatrixUnit();
    const auto foundResolution = _fs->read_as_unsigned<std::int32_t>();

    const auto nValues = readNValues();

    using MT = MatrixType;
    const bool store = chrom1 == chrom2 && (matrix_type == MT::oe || matrix_type == MT::expected) &&
                       foundNorm == wanted_norm && foundUnit == wanted_unit &&
                       foundResolution == wanted_resolution;

    if (store) {
      expectedValues = readExpectedVector(nValues);
      const auto normFactors = readNormalizationFactors(chrom1.id());
      applyNormalizationFactors(expectedValues, normFactors);
    } else {
      discardExpectedVector(nValues);
      discardNormalizationFactors(chrom1.id());
    }
  }
  return expectedValues;
}

[[nodiscard]] inline balancing::Weights default_initialize_weight_vector(
    const Chromosome &chrom, const balancing::Method &norm, std::uint32_t resolution) {
  const auto filler_weight =
      norm == balancing::Method::NONE() ? 1.0 : std::numeric_limits<double>::quiet_NaN();

  const auto num_bins1 = (chrom.size() + resolution - 1) / resolution;
  return {filler_weight, num_bins1, balancing::Weights::Type::DIVISIVE};
}

inline void HiCFileReader::read_footer_norm(const Chromosome &chrom1, const Chromosome &chrom2,
                                            const balancing::Method &wanted_norm,
                                            MatrixUnit wanted_unit, std::uint32_t wanted_resolution,
                                            std::shared_ptr<balancing::Weights> &weights1,
                                            std::shared_ptr<balancing::Weights> &weights2) {
  assert(weights1);
  assert(weights2);
  if (!weights1->empty() && !weights2->empty()) {
    return;
  }

  if (_fs->tellg() == static_cast<std::streampos>(_fs->size())) {
    if (weights1->empty()) {
      *weights1 = default_initialize_weight_vector(chrom1, wanted_norm, wanted_resolution);
    }
    if (weights2->empty()) {
      *weights2 = default_initialize_weight_vector(chrom2, wanted_norm, wanted_resolution);
    }
    return;
  }

  // Index of normalization vectors
  const auto nEntries = _fs->read<std::int32_t>();
  bool norm_found = false;
  for (std::int32_t i = 0; i < nEntries; i++) {
    const auto foundNorm = readNormalizationMethod();
    const auto foundChrom = _fs->read_as_unsigned<std::int32_t>();
    const auto foundUnit = readMatrixUnit();

    const auto foundResolution = _fs->read_as_unsigned<std::int32_t>();
    const auto filePosition = _fs->read<std::int64_t>();
    const auto sizeInBytes = version() > 8 ? _fs->read<std::int64_t>()
                                           : static_cast<std::int64_t>(_fs->read<std::int32_t>());

    norm_found |= foundNorm == wanted_norm && foundUnit == wanted_unit &&
                  foundResolution == wanted_resolution;

    const auto store1 = weights1->empty() && foundChrom == chrom1.id() &&
                        foundNorm == wanted_norm && foundUnit == wanted_unit &&
                        foundResolution == wanted_resolution;
    if (store1) {
      const auto numBins =
          static_cast<std::size_t>((chrom1.size() + wanted_resolution - 1) / wanted_resolution);
      const auto currentPos = static_cast<std::int64_t>(_fs->tellg());
      *weights1 = balancing::Weights{
          readNormalizationVector(indexEntry{filePosition, sizeInBytes}, numBins),
          balancing::Weights::Type::DIVISIVE};
      _fs->seekg(currentPos);
    }

    const auto store2 = weights2->empty() && foundChrom == chrom2.id() &&
                        foundNorm == wanted_norm && foundUnit == wanted_unit &&
                        foundResolution == wanted_resolution;
    if (store2) {
      const auto numBins =
          static_cast<std::size_t>((chrom2.size() + wanted_resolution - 1) / wanted_resolution);
      const auto currentPos = static_cast<std::int64_t>(_fs->tellg());
      *weights2 = balancing::Weights{
          readNormalizationVector(indexEntry{filePosition, sizeInBytes}, numBins),
          balancing::Weights::Type::DIVISIVE};
      _fs->seekg(currentPos);
    }
  }

  if (!norm_found) {
    throw std::runtime_error(fmt::format(FMT_STRING("unable to read \"{}\" weights"), wanted_norm));
  }

  if (weights1->empty()) {
    *weights1 = default_initialize_weight_vector(chrom1, wanted_norm, wanted_resolution);
  }
  if (weights2->empty()) {
    *weights2 = default_initialize_weight_vector(chrom1, wanted_norm, wanted_resolution);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
inline HiCFooter HiCFileReader::read_footer(const Chromosome &chrom1, const Chromosome &chrom2,
                                            MatrixType matrix_type,
                                            const balancing::Method &wanted_norm,
                                            MatrixUnit wanted_unit, std::uint32_t wanted_resolution,
                                            std::shared_ptr<balancing::Weights> &weights1,
                                            std::shared_ptr<balancing::Weights> &weights2) {
  assert(chrom1 <= chrom2);
  assert(std::find(_header->resolutions.begin(), _header->resolutions.end(), wanted_resolution) !=
         _header->resolutions.end());
  assert(!!weights1);
  assert(!!weights2);
  if (chrom1 == chrom2) {
    assert(weights1 == weights2);
  } else {
    assert(weights1 != weights2);
  }

  using MT = MatrixType;
  using NM = balancing::Method;

  // clang-format off
    HiCFooterMetadata metadata{
        _fs->path(),
        matrix_type,
        wanted_norm,
        wanted_unit,
        wanted_resolution,
        chrom1,
        chrom2
    };
  // clang-format on

  auto try_init_weights = [&]() {
    if (weights1->empty()) {
      *weights1 = default_initialize_weight_vector(chrom1, wanted_norm, wanted_resolution);
    }
    if (weights2->empty()) {
      *weights2 = default_initialize_weight_vector(chrom2, wanted_norm, wanted_resolution);
    }
  };

  const auto key = fmt::format(FMT_COMPILE("{}_{}"), chrom1.id(), chrom2.id());

  _fs->seekg(masterOffset());

  metadata.matrixMetadataOffset = read_footer_file_offset(key);
  if (metadata.matrixMetadataOffset == -1) {
    try_init_weights();
    return {Index{}, std::move(metadata), {}, weights1, weights2};
  }

  const auto file_offset = _fs->tellg();
  // NOTE: we read then move index to workaround assertion failures when compiling under MSVC
  auto index = read_index(metadata.matrixMetadataOffset, metadata.chrom1, metadata.chrom2,
                          metadata.unit, metadata.resolution);
  _fs->seekg(static_cast<std::int64_t>(file_offset));

  if ((matrix_type == MT::observed && wanted_norm == NM::NONE()) ||
      ((matrix_type == MT::oe || matrix_type == MT::expected) && wanted_norm == NM::NONE() &&
       chrom1 != chrom2)) {
    // no need to read wanted_norm vector index
    try_init_weights();
    return {std::move(index), std::move(metadata), {}, weights1, weights2};
  }

  auto expectedValues = read_footer_expected_values(chrom1, chrom2, matrix_type, wanted_norm,
                                                    wanted_unit, wanted_resolution);
  if (chrom1 == chrom2 && (matrix_type == MT::oe || matrix_type == MT::expected) &&
      wanted_norm == NM::NONE()) {
    if (expectedValues.empty()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to find expected values for {}:{} at {} ({})"),
                      chrom1.name(), chrom2.name(), wanted_resolution, wanted_unit));
    }
    try_init_weights();
    return {std::move(index), std::move(metadata), std::move(expectedValues), weights1, weights2};
  }

  expectedValues = read_footer_expected_values_norm(chrom1, chrom2, matrix_type, wanted_norm,
                                                    wanted_unit, wanted_resolution);
  if (chrom1 == chrom2 && (matrix_type == MT::oe || matrix_type == MT::expected) &&
      wanted_norm != NM::NONE()) {
    if (expectedValues.empty()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to find normalization factors for {}:{} at "
                                 "{} ({})"),
                      chrom1.name(), chrom2.name(), wanted_resolution, wanted_unit));
    }
  }

  read_footer_norm(chrom1, chrom2, wanted_norm, wanted_unit, wanted_resolution, weights1, weights2);

  return {std::move(index), std::move(metadata), std::move(expectedValues), weights1, weights2};
}

inline std::vector<balancing::Method> HiCFileReader::list_avail_normalizations(
    MatrixType matrix_type, MatrixUnit wanted_unit, std::uint32_t wanted_resolution) {
  if (version() >= 9) {  // NOLINT(*-avoid-magic-numbers)
    return list_avail_normalizations_v9();
  }

  phmap::flat_hash_set<balancing::Method> methods{};
  _fs->seekg(masterOffset());
  [[maybe_unused]] const auto offset = read_footer_file_offset("1_1");
  assert(offset != -1);

  const auto &chrom = _header->chromosomes.longest_chromosome();
  std::ignore = read_footer_expected_values(chrom, chrom, matrix_type, balancing::Method::NONE(),
                                            wanted_unit, wanted_resolution);
  if (_fs->tellg() == static_cast<std::streampos>(_fs->size())) {
    return {};
  }

  std::ignore = read_footer_expected_values_norm(
      chrom, chrom, matrix_type, balancing::Method::NONE(), wanted_unit, wanted_resolution);
  if (_fs->tellg() == static_cast<std::streampos>(_fs->size())) {
    return {};
  }

  const auto nNormVectors = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nNormVectors; i++) {
    const auto foundNorm = readNormalizationMethod();
    methods.emplace(foundNorm);
    [[maybe_unused]] const auto chrIdx = _fs->read<std::int32_t>();
    [[maybe_unused]] const auto foundUnit = readMatrixUnit();
    [[maybe_unused]] const auto foundResolution = _fs->read_as_unsigned<std::int32_t>();
    [[maybe_unused]] const auto position = _fs->read<std::int64_t>();
    [[maybe_unused]] const auto nBytes = _fs->read<std::int32_t>();
  }

  std::vector<balancing::Method> methods_{methods.size()};
  std::copy(methods.begin(), methods.end(), methods_.begin());
  std::sort(methods_.begin(), methods_.end(),
            [&](const auto &m1, const auto &m2) { return m1.to_string() < m2.to_string(); });
  return methods_;
}

inline std::vector<balancing::Method> HiCFileReader::list_avail_normalizations_v9() {
  if (_header->normVectorIndexPosition <= 0) {
    return {};
  }
  phmap::flat_hash_set<balancing::Method> methods{};
  _fs->seekg(_header->normVectorIndexPosition);
  const auto nNormVectors = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nNormVectors; i++) {
    const auto foundNorm = readNormalizationMethod();
    methods.emplace(foundNorm);
    [[maybe_unused]] const auto chrIdx = _fs->read<std::int32_t>();
    [[maybe_unused]] const auto foundUnit = readMatrixUnit();
    [[maybe_unused]] const auto foundResolution = _fs->read_as_unsigned<std::int32_t>();
    [[maybe_unused]] const auto position = _fs->read<std::int64_t>();
    [[maybe_unused]] const auto nBytes = _fs->read<std::int64_t>();
  }

  std::vector<balancing::Method> methods_{methods.size()};
  std::copy(methods.begin(), methods.end(), methods_.begin());
  std::sort(methods_.begin(), methods_.end(),
            [&](const auto &m1, const auto &m2) { return m1.to_string() < m2.to_string(); });
  return methods_;
}

}  // namespace hictk::hic::internal
