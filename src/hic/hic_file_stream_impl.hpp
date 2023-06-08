// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>
#include <utility>

#include "hictk/hic/common.hpp"
#include "hictk/hic/filestream.hpp"

namespace hictk::internal {

inline HiCFileStream::HiCFileStream(std::string url)
    : _fs(std::make_shared<filestream::FileStream>(HiCFileStream::openStream(std::move(url)))),
      _header(std::make_shared<const HiCHeader>(HiCFileStream::readHeader(*_fs))) {}

inline filestream::FileStream HiCFileStream::openStream(std::string url) {
#ifdef HICXX_USE_CURL
  const auto isRemoteFile = StartsWith(url, "http");
  try {
    try {
      if (isRemoteFile) {
        constexpr std::size_t defaulChunkSize = 64 * 1024;
        return filestream::FileStream::Remote(url, defaulChunkSize, "hictk");
      }
      return filestream::FileStream::Local(url);

    } catch (const std::system_error &e) {
      throw std::runtime_error(std::string(e.what()) + ": " + e.code().message());
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Failed to open {} : {}"),
                                         isRemoteFile ? "remote file" : "file", e.what()));
  }
#else
  try {
    return filestream::FileStream(url);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Failed to open file: {}"), e.what()));
  }
#endif
}

inline const std::string &HiCFileStream::url() const noexcept { return _fs->url(); }
inline const HiCHeader &HiCFileStream::header() const noexcept { return *_header; }

inline bool HiCFileStream::isLocal() const noexcept { return !isRemote(); }

inline bool HiCFileStream::isRemote() const noexcept {
#ifdef HICXX_USE_CURL
  return _fs->is_remote();
#else
  return false;
#endif
}

inline std::int32_t HiCFileStream::version() const noexcept {
  assert(_header->version != -1);
  return _header->version;
}

inline void HiCFileStream::discardExpectedVector(std::int64_t nValues) {
  const std::int64_t elementSize = version() > 8 ? sizeof(float) : sizeof(double);
  _fs->seekg(nValues * elementSize, std::ios::cur);
}

inline std::vector<double> HiCFileStream::readExpectedVector(std::int64_t nValues) {
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

inline std::vector<double> HiCFileStream::readNormalizationFactors(std::int32_t wantedChrom) {
  const auto nFactors = _fs->read<std::int32_t>();
  std::vector<double> normFactors{};
  auto readFactor = [this]() {
    if (version() > 8) {
      return static_cast<double>(_fs->read<float>());
    }
    return _fs->read<double>();
  };

  for (auto i = 0; i < nFactors; ++i) {
    const auto foundChrom = _fs->read<std::int32_t>();
    const auto v = readFactor();
    if (foundChrom == wantedChrom) {
      normFactors.push_back(v);
    }
  }
  return normFactors;
}

inline void HiCFileStream::applyNormalizationFactors(std::vector<double> &expectedValues,
                                                     const std::vector<double> &normFactors) {
  if (normFactors.empty() || expectedValues.empty()) {
    return;
  }
  for (const auto factor : normFactors) {
    std::transform(expectedValues.begin(), expectedValues.end(), expectedValues.begin(),
                   [&](auto n) { return n / factor; });
  }
}
inline std::vector<double> HiCFileStream::readNormalizationVector(indexEntry cNormEntry,
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

inline void HiCFileStream::discardNormalizationFactors(std::int32_t wantedChrom) {
  std::ignore = readNormalizationFactors(wantedChrom);
}

inline MatrixType HiCFileStream::readMatrixType(filestream::FileStream &fs, std::string &buff) {
  fs.getline(buff, '\0');
  return ParseMatrixTypeStr(buff);
}

inline NormalizationMethod HiCFileStream::readNormalizationMethod(filestream::FileStream &fs,
                                                                  std::string &buff) {
  fs.getline(buff, '\0');
  return ParseNormStr(buff);
}

inline MatrixUnit HiCFileStream::readMatrixUnit(filestream::FileStream &fs, std::string &buff) {
  fs.getline(buff, '\0');
  return ParseUnitStr(buff);
}

inline MatrixType HiCFileStream::readMatrixType() {
  return HiCFileStream::readMatrixType(*_fs, _strbuff);
}

inline NormalizationMethod HiCFileStream::readNormalizationMethod() {
  return HiCFileStream::readNormalizationMethod(*_fs, _strbuff);
}

inline MatrixUnit HiCFileStream::readMatrixUnit() {
  return HiCFileStream::readMatrixUnit(*_fs, _strbuff);
}

inline std::int64_t HiCFileStream::readNValues() {
  if (version() > 8) {
    return _fs->read<std::int64_t>();
  }
  return _fs->read<std::int32_t>();
}

inline bool HiCFileStream::checkMagicString(filestream::FileStream &fs) {
  return fs.getline('\0') == "HIC";
}

inline std::int64_t HiCFileStream::masterOffset() const noexcept {
  return _header->masterIndexOffset;
}

inline auto HiCFileStream::init_decompressor() -> Decompressor {
  Decompressor zs(libdeflate_alloc_decompressor(),
                  [](auto *ptr) { libdeflate_free_decompressor(ptr); });
  if (!zs) {
    throw std::runtime_error("failed to initialize zlib decompression stream");
  }

  return zs;
}

inline void HiCFileStream::readBlockMap(std::int64_t fileOffset,
                                        [[maybe_unused]] const chromosome &chrom1,
                                        [[maybe_unused]] const chromosome &chrom2,
                                        MatrixUnit wantedUnit, std::int64_t wantedResolution,
                                        BlockMap &buffer) {
  _fs->seekg(fileOffset);
  auto &blockMap = buffer.blocks;
  blockMap.clear();

  [[maybe_unused]] const auto c1i = _fs->read<std::int32_t>();
  [[maybe_unused]] const auto c2i = _fs->read<std::int32_t>();
  const auto numResolutions = _fs->read<std::int32_t>();

  assert(c1i == chrom1.index);
  assert(c2i == chrom2.index);

  for (std::int32_t i = 0; i < numResolutions; ++i) {
    const auto foundUnit = readMatrixUnit();
    std::ignore = _fs->read<std::int32_t>();  // oldIndex
    const auto sumCount = _fs->read<float>();
    std::ignore = _fs->read<float>();  // occupiedCellCount
    std::ignore = _fs->read<float>();  // stdDev
    std::ignore = _fs->read<float>();  // percent95

    const auto foundResolution = static_cast<std::int64_t>(_fs->read<std::int32_t>());
    const auto blockBinCount = _fs->read<std::int32_t>();
    const auto blockColumnCount = _fs->read<std::int32_t>();

    const auto nBlocks = static_cast<std::int64_t>(_fs->read<std::int32_t>());

    if (wantedUnit == foundUnit && wantedResolution == foundResolution) {
      for (std::int64_t j = 0; j < nBlocks; ++j) {
        const auto key = _fs->read<std::int32_t>();
        indexEntry index{_fs->read<std::int64_t>(), _fs->read<std::int32_t>()};
        assert(index.position + index.size < static_cast<std::int64_t>(_fs->size()));
        blockMap.emplace(key, std::move(index));
      }
      buffer.blockBinCount = blockBinCount;
      buffer.blockColumnCount = blockColumnCount;
      buffer.sumCount = static_cast<double>(sumCount);
      return;
    }

    constexpr std::int64_t blockSize = sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t);
    _fs->seekg(nBlocks * blockSize, std::ios::cur);
  }

  throw std::runtime_error(
      fmt::format(FMT_STRING("Unable to find block map for unit {} and resolution {}"), wantedUnit,
                  wantedResolution));
}

inline bool HiCFileStream::checkMagicString() { return checkMagicString(*_fs); }

// reads the header, storing the positions of the normalization vectors and returning the
// masterIndexPosition pointer
inline HiCHeader HiCFileStream::readHeader(filestream::FileStream &fs) {
  if (!checkMagicString(fs)) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("Hi-C magic string is missing. {} does not appear to be a hic file"), fs.url()));
  }

  HiCHeader header{fs.url()};

  fs.read(header.version);
  if (header.version < 8) {
    throw std::runtime_error(
        fmt::format(FMT_STRING(".hic version 7 and older ar no longer supported. Found version {}"),
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

  // reading and ignoring attribute-value dictionary
  for (std::int32_t i = 0; i < nAttributes; i++) {
    std::ignore = fs.getline('\0');  // key
    std::ignore = fs.getline('\0');  // value
  }

  // Read chromosomes
  const auto numChromosomes = fs.read<std::int32_t>();
  std::generate_n(std::inserter(header.chromosomes, header.chromosomes.begin()), numChromosomes,
                  [&]() {
                    chromosome chrom{};
                    chrom.index = static_cast<std::int32_t>(header.chromosomes.size());
                    fs.getline(chrom.name, '\0');
                    if (header.version > 8) {
                      fs.read(chrom.length);
                    } else {
                      chrom.length = static_cast<std::int64_t>(fs.read<std::int32_t>());
                    }

                    return std::make_pair(chrom.name, chrom);
                  });

  if (header.chromosomes.empty()) {
    throw std::runtime_error("unable to read chromosomes");
  }

  // Read resolutions
  const auto numResolutions = static_cast<std::size_t>(fs.read<std::int32_t>());
  if (numResolutions == 0) {
    throw std::runtime_error("unable to read the list of available resolutions");
  }
  header.resolutions.resize(numResolutions);
  std::generate(header.resolutions.begin(), header.resolutions.end(), [&]() {
    const auto res = fs.read<std::int32_t>();
    assert(res > 0);
    return res;
  });

  return header;
}

inline void HiCFileStream::readAndInflate(indexEntry idx, std::string &plainTextBuffer) {
  try {
    // _strbuff is used to store compressed data
    // plainTextBuffer is used to store decompressed data
    assert(_decompressor);
    assert(idx.size > 0);
    const auto buffSize = static_cast<std::size_t>(idx.size);

    plainTextBuffer.reserve(buffSize * 3);
    plainTextBuffer.resize(plainTextBuffer.capacity());

    _fs->seekg(idx.position);
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
                                         idx.position, e.what()));
  }
}

inline bool HiCFileStream::checkMagicString(std::string url) noexcept {
  try {
    filestream::FileStream fs(HiCFileStream::openStream(std::move(url)));
    return HiCFileStream::checkMagicString(fs);
  } catch (...) {
    return false;
  }
}

inline HiCFooter HiCFileStream::readFooter(const std::int32_t chrom1_id,
                                           const std::int32_t chrom2_id,
                                           const MatrixType wantedMatrixType,
                                           const NormalizationMethod wantedNorm,
                                           const MatrixUnit wantedUnit,
                                           const std::int32_t wantedResolution) {
  assert(chrom1_id <= chrom2_id);
  assert(std::find(_header->resolutions.begin(), _header->resolutions.end(), wantedResolution) !=
         _header->resolutions.end());

  using MT = MatrixType;
  using NM = NormalizationMethod;

  // clang-format off
    HiCFooter footer{
        HiCFooterMetadata{_fs->url(),
                          wantedMatrixType,
                          wantedNorm,
                          wantedUnit,
                          wantedResolution,
                          _header->getChromosome(chrom1_id),
                          _header->getChromosome(chrom2_id)}
        };
  // clang-format on

  auto &metadata = footer.metadata();
  auto &expectedValues = footer.expectedValues();
  auto &c1Norm = footer.c1Norm();
  auto &c2Norm = footer.c2Norm();

  const auto key = std::to_string(chrom1_id) + "_" + std::to_string(chrom2_id);

  _fs->seekg(masterOffset());
  std::ignore = readNValues();  // nBytes

  auto nEntries = _fs->read<std::int32_t>();
  for (int i = 0; i < nEntries; i++) {
    const auto strbuff = _fs->getline('\0');
    const auto fpos = _fs->read<std::int64_t>();
    std::ignore = _fs->read<std::int32_t>();  // sizeInBytes
    if (strbuff == key) {
      metadata.fileOffset = fpos;
    }
  }
  if (metadata.fileOffset == -1) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("unable to find interactions for {}:{} at {} ({}): unable to read file offset"),
        _header->getChromosome(chrom1_id).name, _header->getChromosome(chrom2_id).name,
        wantedResolution, wantedUnit));
  }

  if ((wantedMatrixType == MT::observed && wantedNorm == NM::NONE) ||
      ((wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) && wantedNorm == NM::NONE &&
       chrom1_id != chrom2_id)) {
    return footer;  // no need to read wantedNorm vector index
  }

  // read in and ignore expected value maps; don't store; reading these to
  // get to wantedNorm vector index
  auto nExpectedValues = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nExpectedValues; i++) {
    const auto foundUnit = readMatrixUnit();
    const auto foundResolution = _fs->read<std::int32_t>();
    const auto nValues = readNValues();

    bool store = chrom1_id == chrom2_id &&
                 (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
                 wantedNorm == NM::NONE && foundUnit == wantedUnit &&
                 foundResolution == wantedResolution;

    if (store) {
      expectedValues = readExpectedVector(nValues);
      const auto normFactors = readNormalizationFactors(chrom1_id);
      applyNormalizationFactors(expectedValues, normFactors);

    } else {
      discardExpectedVector(nValues);
      discardNormalizationFactors(chrom1_id);
    }
  }

  if (chrom1_id == chrom2_id && (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
      wantedNorm == NM::NONE) {
    if (expectedValues.empty()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to find expected values for {}:{} at {} ({})"),
                      _header->getChromosome(chrom1_id).name,
                      _header->getChromosome(chrom2_id).name, wantedResolution, wantedUnit));
    }
    return footer;
  }

  nExpectedValues = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nExpectedValues; i++) {
    const auto foundNorm = readNormalizationMethod();
    const auto foundUnit = readMatrixUnit();
    const auto foundResolution = _fs->read<std::int32_t>();

    const auto nValues = readNValues();
    bool store = chrom1_id == chrom2_id &&
                 (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
                 foundNorm == wantedNorm && foundUnit == wantedUnit &&
                 foundResolution == wantedResolution;

    if (store) {
      expectedValues = readExpectedVector(nValues);
      const auto normFactors = readNormalizationFactors(chrom1_id);
      applyNormalizationFactors(expectedValues, normFactors);
    } else {
      discardExpectedVector(nValues);
      discardNormalizationFactors(chrom1_id);
    }
  }

  if (chrom1_id == chrom2_id && (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
      wantedNorm != NM::NONE) {
    if (expectedValues.empty()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("unable to find expected values normalization factors for {}:{} at "
                     "{} ({})"),
          _header->getChromosome(chrom1_id).name, _header->getChromosome(chrom2_id).name,
          wantedResolution, wantedUnit));
    }
  }

  // Index of normalization vectors
  nEntries = _fs->read<std::int32_t>();
  for (std::int32_t i = 0; i < nEntries; i++) {
    const auto foundNorm = readNormalizationMethod();
    const auto foundChrom = _fs->read<std::int32_t>();
    const auto foundUnit = readMatrixUnit();

    const auto foundResolution = _fs->read<std::int32_t>();
    const auto filePosition = _fs->read<std::int64_t>();
    const auto sizeInBytes = version() > 8 ? _fs->read<std::int64_t>()
                                           : static_cast<std::int64_t>(_fs->read<std::int32_t>());
    if (foundChrom == chrom1_id && foundNorm == wantedNorm && foundUnit == wantedUnit &&
        foundResolution == wantedResolution) {
      const auto numBins = static_cast<std::size_t>(
          (footer.chrom1().length + wantedResolution - 1) / wantedResolution);
      const auto currentPos = static_cast<std::int64_t>(this->_fs->tellg());
      c1Norm = readNormalizationVector(indexEntry{filePosition, sizeInBytes}, numBins);
      _fs->seekg(currentPos);
    }
    if (chrom1_id != chrom2_id && foundChrom == chrom2_id && foundNorm == wantedNorm &&
        foundUnit == wantedUnit && foundResolution == wantedResolution) {
      const auto numBins = static_cast<std::size_t>(
          (footer.chrom2().length + wantedResolution - 1) / wantedResolution);
      const auto currentPos = static_cast<std::int64_t>(this->_fs->tellg());
      c2Norm = readNormalizationVector(indexEntry{filePosition, sizeInBytes}, numBins);
      _fs->seekg(currentPos);
    }
  }

  if (footer.c1Norm().empty() && footer.c2Norm().empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find {} normalization vectors for {}:{} at {} ({})"),
                    wantedNorm, _header->getChromosome(chrom1_id).name,
                    _header->getChromosome(chrom2_id).name, wantedResolution, wantedUnit));
  }

  if (footer.c1Norm().empty() || footer.c2Norm().empty()) {
    const auto chrom_id = footer.c1Norm().empty() ? chrom1_id : chrom2_id;
    throw std::runtime_error(fmt::format(
        FMT_STRING("unable to find {} normalization vector for {} at {} ({})"), wantedNorm,
        _header->getChromosome(chrom_id).name, wantedResolution, wantedUnit));
  }

  return footer;
}

}  // namespace hictk::internal
