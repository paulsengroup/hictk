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
#include "hictk/hic/common.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/reference.hpp"

namespace hictk::hic::internal {

inline HiCFileReader::HiCFileReader(std::string url)
    : _fs(std::make_shared<filestream::FileStream>(HiCFileReader::openStream(std::move(url)))),
      _header(std::make_shared<const HiCHeader>(HiCFileReader::readHeader(*_fs))) {}

inline filestream::FileStream HiCFileReader::openStream(std::string url) {
  try {
    return filestream::FileStream(url);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Failed to open file: {}"), e.what()));
  }
}

inline const std::string &HiCFileReader::url() const noexcept { return _fs->url(); }
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
  if (version() > 8) {
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
    if (version() > 8) {
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
  if (version() > 8) {
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

inline MatrixType HiCFileReader::readMatrixType(filestream::FileStream &fs, std::string &buff) {
  fs.getline(buff, '\0');
  return ParseMatrixTypeStr(buff);
}

inline balancing::Method HiCFileReader::readNormalizationMethod(filestream::FileStream &fs,
                                                                std::string &buff) {
  fs.getline(buff, '\0');
  return balancing::Method{buff};
}

inline MatrixUnit HiCFileReader::readMatrixUnit(filestream::FileStream &fs, std::string &buff) {
  fs.getline(buff, '\0');
  return ParseUnitStr(buff);
}

inline MatrixType HiCFileReader::readMatrixType() {
  return HiCFileReader::readMatrixType(*_fs, _strbuff);
}

inline balancing::Method HiCFileReader::readNormalizationMethod() {
  return HiCFileReader::readNormalizationMethod(*_fs, _strbuff);
}

inline MatrixUnit HiCFileReader::readMatrixUnit() {
  return HiCFileReader::readMatrixUnit(*_fs, _strbuff);
}

inline std::int64_t HiCFileReader::readNValues() {
  if (version() > 8) {
    return _fs->read<std::int64_t>();
  }
  return _fs->read<std::int32_t>();
}

inline bool HiCFileReader::checkMagicString(filestream::FileStream &fs) {
  return fs.getline('\0') == "HIC";
}

inline std::int64_t HiCFileReader::masterOffset() const noexcept {
  return _header->masterIndexOffset;
}

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

    Index::BlkIdxBuffer buffer;
    buffer.reserve(nBlocks);
    if (wantedUnit == foundUnit && wantedResolution == foundResolution) {
      for (std::size_t j = 0; j < nBlocks; ++j) {
        const auto block_id = static_cast<std::size_t>(_fs->read<std::int32_t>());
        const auto position = static_cast<std::size_t>(_fs->read<std::int64_t>());
        const auto size = static_cast<std::size_t>(_fs->read<std::int32_t>());
        assert(position + size < _fs->size());
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

// reads the header, storing the positions of the normalization vectors and returning the
// masterIndexPosition pointer
inline HiCHeader HiCFileReader::readHeader(filestream::FileStream &fs) {
  if (!checkMagicString(fs)) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Hi-C magic string is missing. {} does not appear to be a hic file"), fs.url()));
  }

  HiCHeader header{fs.url()};

  fs.read(header.version);
  if (header.version < 6) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(".hic version 5 and older are no longer supported. Found version {}"),
        header.version));
  }
  fs.read(header.masterIndexOffset);
  if (header.masterIndexOffset < 0 ||
      header.masterIndexOffset >= static_cast<std::int64_t>(fs.size())) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("file appears to be corrupted: expected master index offset to "
                               "be between 0 and {}, found {}"),
                    fs.size(), header.masterIndexOffset));
  }

  fs.getline(header.genomeID, '\0');
  if (header.genomeID.empty()) {
    header.genomeID = "unknown";
  }

  if (header.version > 8) {
    fs.read(header.nviPosition);
    fs.read(header.nviLength);
  }

  const auto nAttributes = fs.read<std::int32_t>();

  // reading attribute-value dictionary
  for (std::int32_t i = 0; i < nAttributes; i++) {
    auto key = fs.getline('\0');    // key
    auto value = fs.getline('\0');  // value
    header.attributes.emplace(std::move(key), std::move(value));
  }

  // Read chromosomes
  auto numChromosomes = static_cast<std::uint32_t>(fs.read<std::int32_t>());
  std::vector<std::string> chrom_names(numChromosomes);
  std::vector<std::uint32_t> chrom_sizes(numChromosomes);
  for (std::size_t i = 0; i < chrom_names.size(); ++i) {
    fs.getline(chrom_names[i], '\0');
    chrom_sizes[i] = static_cast<std::uint32_t>(
        header.version > 8 ? fs.read<std::int64_t>()
                           : static_cast<std::int64_t>(fs.read<std::int32_t>()));
  }

  if (chrom_names.empty()) {
    throw std::runtime_error("unable to read chromosomes");
  }

  header.chromosomes = Reference(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());

  // Read resolutions
  const auto numResolutions = static_cast<std::size_t>(fs.read<std::int32_t>());
  if (numResolutions == 0) {
    throw std::runtime_error("unable to read the list of available resolutions");
  }
  header.resolutions.resize(numResolutions);
  std::generate(header.resolutions.begin(), header.resolutions.end(), [&]() {
    const auto res = fs.read<std::int32_t>();
    assert(res > 0);
    return static_cast<std::uint32_t>(res);
  });

  return header;
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
    filestream::FileStream fs(HiCFileReader::openStream(std::move(url)));
    return HiCFileReader::checkMagicString(fs);
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
    const auto fpos = _fs->read<std::int64_t>();
    std::ignore = _fs->read<std::int32_t>();  // sizeInBytes
    if (strbuff == key) {
      pos = fpos;
    }
  }
  return pos;
}

inline std::vector<double> HiCFileReader::read_footer_expected_values(
    std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
    balancing::Method wanted_norm, MatrixUnit wanted_unit, std::uint32_t wanted_resolution) {
  std::vector<double> expectedValues{};
  auto nExpectedValues = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nExpectedValues; ++i) {
    const auto foundUnit = readMatrixUnit();
    const auto foundResolution = _fs->read_as_unsigned<std::int32_t>();
    const auto nValues = readNValues();

    using MT = MatrixType;
    using NM = balancing::Method;
    bool store = chrom1_id == chrom2_id && (matrix_type == MT::oe || matrix_type == MT::expected) &&
                 wanted_norm == NM::NONE() && foundUnit == wanted_unit &&
                 foundResolution == wanted_resolution;

    if (store) {
      expectedValues = readExpectedVector(nValues);
      const auto normFactors = readNormalizationFactors(chrom1_id);
      applyNormalizationFactors(expectedValues, normFactors);

    } else {
      discardExpectedVector(nValues);
      discardNormalizationFactors(chrom1_id);
    }
  }
  return expectedValues;
}

inline std::vector<double> HiCFileReader::read_footer_expected_values_norm(
    std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
    balancing::Method wanted_norm, MatrixUnit wanted_unit, std::uint32_t wanted_resolution) {
  if (_fs->tellg() == _fs->size()) {
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
    bool store = chrom1_id == chrom2_id && (matrix_type == MT::oe || matrix_type == MT::expected) &&
                 foundNorm == wanted_norm && foundUnit == wanted_unit &&
                 foundResolution == wanted_resolution;

    if (store) {
      expectedValues = readExpectedVector(nValues);
      const auto normFactors = readNormalizationFactors(chrom1_id);
      applyNormalizationFactors(expectedValues, normFactors);
    } else {
      discardExpectedVector(nValues);
      discardNormalizationFactors(chrom1_id);
    }
  }
  return expectedValues;
}

inline void HiCFileReader::read_footer_norm(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                            balancing::Method wanted_norm, MatrixUnit wanted_unit,
                                            std::uint32_t wanted_resolution,
                                            const Chromosome &chrom1, const Chromosome &chrom2,
                                            std::shared_ptr<balancing::Weights> &weights1,
                                            std::shared_ptr<balancing::Weights> &weights2) {
  if (_fs->tellg() == _fs->size()) {
    return;
  }

  // Index of normalization vectors
  const auto nEntries = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nEntries; i++) {
    const auto foundNorm = readNormalizationMethod();
    const auto foundChrom = _fs->read_as_unsigned<std::int32_t>();
    const auto foundUnit = readMatrixUnit();

    const auto foundResolution = _fs->read_as_unsigned<std::int32_t>();
    const auto filePosition = _fs->read<std::int64_t>();
    const auto sizeInBytes = version() > 8 ? _fs->read<std::int64_t>()
                                           : static_cast<std::int64_t>(_fs->read<std::int32_t>());

    const auto store1 = !*weights1 && foundChrom == chrom1_id && foundNorm == wanted_norm &&
                        foundUnit == wanted_unit && foundResolution == wanted_resolution;
    if (store1) {
      const auto numBins =
          static_cast<std::size_t>((chrom1.size() + wanted_resolution - 1) / wanted_resolution);
      const auto currentPos = static_cast<std::int64_t>(_fs->tellg());
      *weights1 = balancing::Weights{
          readNormalizationVector(indexEntry{filePosition, sizeInBytes}, numBins),
          balancing::Weights::Type::DIVISIVE};
      _fs->seekg(currentPos);
    }

    const auto store2 = !*weights2 && foundChrom == chrom2_id && foundNorm == wanted_norm &&
                        foundUnit == wanted_unit && foundResolution == wanted_resolution;
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
}

inline HiCFooter HiCFileReader::read_footer(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                            MatrixType matrix_type, balancing::Method wanted_norm,
                                            MatrixUnit wanted_unit, std::uint32_t wanted_resolution,
                                            std::shared_ptr<balancing::Weights> weights1,
                                            std::shared_ptr<balancing::Weights> weights2) {
  assert(chrom1_id <= chrom2_id);
  assert(std::find(_header->resolutions.begin(), _header->resolutions.end(), wanted_resolution) !=
         _header->resolutions.end());

  using MT = MatrixType;
  using NM = balancing::Method;

  // clang-format off
    HiCFooterMetadata metadata{
        _fs->url(),
        matrix_type,
        wanted_norm,
        wanted_unit,
        wanted_resolution,
        _header->chromosomes.at(chrom1_id),
        _header->chromosomes.at(chrom2_id)
    };
  // clang-format on

  const auto key = fmt::format(FMT_COMPILE("{}_{}"), chrom1_id, chrom2_id);

  _fs->seekg(masterOffset());

  metadata.fileOffset = read_footer_file_offset(key);
  if (metadata.fileOffset == -1) {
    return {Index{}, std::move(metadata), {}, std::move(weights1), std::move(weights2)};
  }

  const auto file_offset = _fs->tellg();
  // NOTE: we read then move index to workaround assertion failures when compiling under MSVC
  auto index = read_index(metadata.fileOffset, metadata.chrom1, metadata.chrom2, metadata.unit,
                          metadata.resolution);
  _fs->seekg(static_cast<std::int64_t>(file_offset));

  if ((matrix_type == MT::observed && wanted_norm == NM::NONE()) ||
      ((matrix_type == MT::oe || matrix_type == MT::expected) && wanted_norm == NM::NONE() &&
       chrom1_id != chrom2_id)) {
    // no need to read wanted_norm vector index
    return {std::move(index), std::move(metadata), {}, std::move(weights1), std::move(weights2)};
  }

  auto expectedValues = read_footer_expected_values(chrom1_id, chrom2_id, matrix_type, wanted_norm,
                                                    wanted_unit, wanted_resolution);
  if (chrom1_id == chrom2_id && (matrix_type == MT::oe || matrix_type == MT::expected) &&
      wanted_norm == NM::NONE()) {
    if (expectedValues.empty()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to find expected values for {}:{} at {} ({})"),
                      _header->chromosomes.at(chrom1_id).name(),
                      _header->chromosomes.at(chrom2_id).name(), wanted_resolution, wanted_unit));
    }
    return {std::move(index), std::move(metadata), std::move(expectedValues), std::move(weights1),
            std::move(weights2)};
  }

  expectedValues = read_footer_expected_values_norm(chrom1_id, chrom2_id, matrix_type, wanted_norm,
                                                    wanted_unit, wanted_resolution);
  if (chrom1_id == chrom2_id && (matrix_type == MT::oe || matrix_type == MT::expected) &&
      wanted_norm != NM::NONE()) {
    if (expectedValues.empty()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to find normalization factors for {}:{} at "
                                 "{} ({})"),
                      _header->chromosomes.at(chrom1_id).name(),
                      _header->chromosomes.at(chrom2_id).name(), wanted_resolution, wanted_unit));
    }
  }

  read_footer_norm(chrom1_id, chrom2_id, wanted_norm, wanted_unit, wanted_resolution,
                   metadata.chrom1, metadata.chrom2, weights1, weights2);

  if (!*weights1 && !*weights2) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find {} normalization vectors for {}:{} at {} ({})"),
                    wanted_norm, _header->chromosomes.at(chrom1_id).name(),
                    _header->chromosomes.at(chrom2_id).name(), wanted_resolution, wanted_unit));
  }

  if (!*weights1 || !*weights2) {
    const auto chrom_id = !*weights1 ? chrom1_id : chrom2_id;
    throw std::runtime_error(fmt::format(
        FMT_STRING("unable to find {} normalization vector for {} at {} ({})"), wanted_norm,
        _header->chromosomes.at(chrom_id).name(), wanted_resolution, wanted_unit));
  }

  return {std::move(index), std::move(metadata), std::move(expectedValues), std::move(weights1),
          std::move(weights2)};
}

inline std::vector<balancing::Method> HiCFileReader::list_avail_normalizations(
    MatrixType matrix_type, MatrixUnit wanted_unit, std::uint32_t wanted_resolution) {
  phmap::flat_hash_set<balancing::Method> methods{};
  _fs->seekg(masterOffset());
  [[maybe_unused]] const auto offset = read_footer_file_offset("1_1");
  assert(offset != -1);

  std::ignore = read_footer_expected_values(1, 1, matrix_type, balancing::Method::NONE(),
                                            wanted_unit, wanted_resolution);
  if (_fs->tellg() == _fs->size()) {
    return {};
  }

  const auto nExpectedValues = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nExpectedValues; i++) {
    const auto foundNorm = readNormalizationMethod();
    methods.emplace(foundNorm);
    [[maybe_unused]] const auto foundUnit = readMatrixUnit();
    [[maybe_unused]] const auto foundResolution = _fs->read_as_unsigned<std::int32_t>();

    [[maybe_unused]] const auto nValues = readNValues();

    discardExpectedVector(nValues);
    discardNormalizationFactors(1);
  }

  std::vector<balancing::Method> methods_{methods.size()};
  std::copy(methods.begin(), methods.end(), methods_.begin());
  std::sort(methods_.begin(), methods_.end(),
            [&](const auto &m1, const auto &m2) { return m1.to_string() < m2.to_string(); });
  return methods_;
}

}  // namespace hictk::hic::internal
