// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/file_writer_data_structures.hpp"

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/binary_buffer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

std::string MatrixMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(chr1Idx);
    buffer.write(chr2Idx);
    buffer.write(nResolutions);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a MatrixMetadata object: {}"), e.what()));
  }
  return buffer.get();
}

std::string MatrixBlockMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(blockNumber);
    buffer.write(blockPosition);
    buffer.write(blockSizeBytes);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a MatrixBlockMetadata object: {}"),
        e.what()));
  }

  return buffer.get();
}

bool MatrixBlockMetadata::operator<(const MatrixBlockMetadata &other) const noexcept {
  return blockNumber < other.blockNumber;
}

bool MatrixResolutionMetadata::operator<(const MatrixResolutionMetadata &other) const noexcept {
  if (unit != other.unit) {
    return unit < other.unit;
  }
  return binSize < other.binSize;
}

std::string MatrixResolutionMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(unit);
    buffer.write(resIdx);
    buffer.write(sumCounts);
    buffer.write(occupiedCellCount);
    buffer.write(percent5);
    buffer.write(percent95);
    buffer.write(binSize);
    buffer.write(blockSize);
    buffer.write(blockColumnCount);
    buffer.write(blockCount);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a MatrixResolutionMetadata object: {}"),
        e.what()));
  }

  for (const auto &blk : _block_metadata) {
    std::ignore = blk.serialize(buffer, false);
  }

  return buffer.get();
}

std::string MatrixBodyMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  try {
    std::ignore = matrixMetadata.serialize(buffer, clear);
    for (const auto &metadata : resolutionMetadata) {
      std::ignore = metadata.serialize(buffer, false);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a MatrixBodyMetadata object: {}"),
        e.what()));
  }

  return buffer.get();
}

std::string FooterMasterIndex::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(key);
    buffer.write(position);
    buffer.write(size);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a FooterMasterIndex object: {}"),
        e.what()));
  }

  return buffer.get();
}

std::int64_t ExpectedValuesBlock::nValues() const noexcept {
  return static_cast<std::int64_t>(value.size());
}

std::int32_t ExpectedValuesBlock::nChrScaleFactors() const noexcept {
  assert(chrIndex.size() == chrScaleFactor.size());
  return static_cast<std::int32_t>(chrIndex.size());
}

ExpectedValuesBlock::ExpectedValuesBlock(std::string_view unit_, std::uint32_t bin_size,
                                         const std::vector<double> &weights,
                                         const std::vector<std::uint32_t> &chrom_ids,
                                         const std::vector<double> &scale_factors)
    : unit(std::string{unit_}),
      binSize(static_cast<std::int32_t>(bin_size)),
      value(weights.size()),
      chrIndex(chrom_ids.size()),
      chrScaleFactor(chrom_ids.size()) {
  std::transform(weights.begin(), weights.end(), value.begin(),
                 [](const auto n) { return static_cast<float>(n); });
  std::transform(chrom_ids.begin(), chrom_ids.end(), chrIndex.begin(),
                 [](const auto n) { return static_cast<std::int32_t>(n); });
  std::transform(scale_factors.begin(), scale_factors.end(), chrScaleFactor.begin(),
                 [](const auto n) { return static_cast<float>(n); });
}

bool ExpectedValuesBlock::operator<(const ExpectedValuesBlock &other) const noexcept {
  if (unit != other.unit) {
    return unit < other.unit;
  }

  return binSize < other.binSize;
}

std::string ExpectedValuesBlock::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(unit);
    buffer.write(binSize);
    buffer.write(nValues());
    buffer.write(value);
    buffer.write(nChrScaleFactors());

    assert(chrIndex.size() == chrScaleFactor.size());
    for (std::size_t i = 0; i < chrIndex.size(); ++i) {
      buffer.write(chrIndex[i]);
      buffer.write(chrScaleFactor[i]);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing an ExpectedValuesBlock object: {}"),
        e.what()));
  }

  return buffer.get();
}

ExpectedValuesBlock ExpectedValuesBlock::deserialize(std::streampos offset,
                                                     filestream::FileStream &fs) {
  [[maybe_unused]] const auto lck = fs.lock();
  return unsafe_deserialize(offset, fs);
}

ExpectedValuesBlock ExpectedValuesBlock::unsafe_deserialize(std::streampos offset,
                                                            filestream::FileStream &fs) {
  assert(offset >= 0);
  ExpectedValuesBlock evb{};

  try {
    fs.unsafe_seekg(offset);
    fs.unsafe_getline(evb.unit, '\0');
    fs.unsafe_read(evb.binSize);
    const auto nValues = static_cast<std::size_t>(fs.unsafe_read<std::int64_t>());
    evb.value.resize(nValues);
    fs.unsafe_read(evb.value);
    const auto nChrScaleFactors = static_cast<std::size_t>(fs.unsafe_read<std::int32_t>());
    evb.chrIndex.resize(nChrScaleFactors);
    evb.chrScaleFactor.resize(nChrScaleFactors);

    for (std::size_t i = 0; i < nChrScaleFactors; ++i) {
      evb.chrIndex.emplace_back(fs.unsafe_read<std::int32_t>());
      evb.chrScaleFactor.emplace_back(fs.unsafe_read<float>());
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while deserializing an ExpectedValuesBlock object: {}"),
        e.what()));
  }

  return evb;
}

std::int32_t ExpectedValues::nExpectedValueVectors() const noexcept {
  return static_cast<std::int32_t>(expectedValues().size());
}

const phmap::btree_set<ExpectedValuesBlock> &ExpectedValues::expectedValues() const noexcept {
  return _expected_values;
}

void ExpectedValues::emplace(const ExpectedValuesBlock &evb, bool force_overwrite) {
  auto [it, inserted] = _expected_values.emplace(evb);
  if (!inserted) {
    if (force_overwrite) {
      *it = evb;
    } else {
      throw std::runtime_error(
          fmt::format(FMT_STRING("ExpectedValues already contains vector for {} resolution ({})"),
                      it->binSize, it->unit));
    }
  }
}

std::string ExpectedValues::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(nExpectedValueVectors());

    if (nExpectedValueVectors() == 0) {
      return buffer.get();
    }

    for (const auto &ev : expectedValues()) {
      std::ignore = ev.serialize(buffer, false);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing an ExpectedValues object: {}"), e.what()));
  }

  return buffer.get();
}

ExpectedValues ExpectedValues::deserialize(std::streampos offset, filestream::FileStream &fs) {
  [[maybe_unused]] const auto lck = fs.lock();
  return unsafe_deserialize(offset, fs);
}

ExpectedValues ExpectedValues::unsafe_deserialize(std::streampos offset,
                                                  filestream::FileStream &fs) {
  assert(offset >= 0);
  ExpectedValues evs{};

  try {
    fs.unsafe_seekg(offset);
    const auto nExpectedValueVectors = static_cast<std::size_t>(fs.unsafe_read<std::int32_t>());
    for (std::size_t i = 0; i < nExpectedValueVectors; ++i) {
      evs.emplace(ExpectedValuesBlock::unsafe_deserialize(fs.unsafe_tellg(), fs), true);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while deserializing an ExpectedValues object: {}"),
        e.what()));
  }
  return evs;
}

std::int64_t NormalizedExpectedValuesBlock::nValues() const noexcept {
  return static_cast<std::int64_t>(value.size());
}

std::int32_t NormalizedExpectedValuesBlock::nChrScaleFactors() const noexcept {
  assert(chrIndex.size() == chrScaleFactor.size());
  return static_cast<std::int32_t>(chrIndex.size());
}

NormalizedExpectedValuesBlock::NormalizedExpectedValuesBlock(
    std::string_view type_, std::string_view unit_, std::uint32_t bin_size,
    const std::vector<double> &weights, const std::vector<std::uint32_t> &chrom_ids,
    const std::vector<double> &scale_factors)
    : type(std::string{type_}),
      unit(std::string{unit_}),
      binSize(static_cast<std::int32_t>(bin_size)),
      value(weights.size()),
      chrIndex(chrom_ids.size()),
      chrScaleFactor(chrom_ids.size()) {
  std::transform(weights.begin(), weights.end(), value.begin(),
                 [](const auto n) { return static_cast<float>(n); });
  std::transform(chrom_ids.begin(), chrom_ids.end(), chrIndex.begin(),
                 [](const auto n) { return static_cast<std::int32_t>(n); });
  std::transform(scale_factors.begin(), scale_factors.end(), chrScaleFactor.begin(),
                 [](const auto n) { return static_cast<float>(n); });
}

bool NormalizedExpectedValuesBlock::operator<(
    const NormalizedExpectedValuesBlock &other) const noexcept {
  if (type != other.type) {
    return type < other.type;
  }
  if (unit != other.unit) {
    return unit < other.unit;
  }
  return binSize < other.binSize;
}

std::string NormalizedExpectedValuesBlock::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(type);
    buffer.write(unit);
    buffer.write(binSize);
    buffer.write(nValues());
    buffer.write(value);
    buffer.write(nChrScaleFactors());

    assert(chrIndex.size() == chrScaleFactor.size());
    for (std::size_t i = 0; i < chrIndex.size(); ++i) {
      buffer.write(chrIndex[i]);
      buffer.write(chrScaleFactor[i]);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while serializing a NormalizedExpectedValuesBlock object: {}"),
        e.what()));
  }

  return buffer.get();
}
NormalizedExpectedValuesBlock NormalizedExpectedValuesBlock::deserialize(
    std::streampos offset, filestream::FileStream &fs) {
  [[maybe_unused]] const auto lck = fs.lock();
  return unsafe_deserialize(offset, fs);
}

NormalizedExpectedValuesBlock NormalizedExpectedValuesBlock::unsafe_deserialize(
    std::streampos offset, filestream::FileStream &fs) {
  assert(offset >= 0);
  NormalizedExpectedValuesBlock nevb{};

  try {
    fs.unsafe_seekg(offset);
    fs.unsafe_getline(nevb.type, '\0');
    fs.unsafe_getline(nevb.unit, '\0');
    fs.unsafe_read(nevb.binSize);
    const auto nValues = static_cast<std::size_t>(fs.unsafe_read<std::int64_t>());
    nevb.value.resize(nValues);
    fs.unsafe_read(nevb.value);
    const auto nChrScaleFactors = static_cast<std::size_t>(fs.unsafe_read<std::int32_t>());
    nevb.chrIndex.resize(nChrScaleFactors);
    nevb.chrScaleFactor.resize(nChrScaleFactors);

    for (std::size_t i = 0; i < nChrScaleFactors; ++i) {
      nevb.chrIndex.emplace_back(fs.unsafe_read<std::int32_t>());
      nevb.chrScaleFactor.emplace_back(fs.unsafe_read<float>());
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while deserializing a NormalizedExpectedValuesBlock object: {}"),
        e.what()));
  }

  return nevb;
}

std::int32_t NormalizedExpectedValues::nNormExpectedValueVectors() const noexcept {
  return static_cast<std::int32_t>(_normalized_expected_values.size());
}

const phmap::btree_set<NormalizedExpectedValuesBlock> &
NormalizedExpectedValues::normExpectedValues() const noexcept {
  return _normalized_expected_values;
}

void NormalizedExpectedValues::emplace(const NormalizedExpectedValuesBlock &evb,
                                       bool force_overwrite) {
  auto [it, inserted] = _normalized_expected_values.emplace(evb);
  if (!inserted) {
    if (force_overwrite) {
      *it = evb;
    } else {
      throw std::runtime_error(fmt::format(
          FMT_STRING("NormalizedExpectedValues already contains {} vector for {} resolution ({})"),
          it->type, it->binSize, it->unit));
    }
  }
}

std::string NormalizedExpectedValues::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(nNormExpectedValueVectors());
    for (const auto &nev : _normalized_expected_values) {
      std::ignore = nev.serialize(buffer, false);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a NormalizedExpectedValues object: {}"),
        e.what()));
  }

  return buffer.get();
}

NormalizedExpectedValues NormalizedExpectedValues::deserialize(std::streampos offset,
                                                               filestream::FileStream &fs) {
  [[maybe_unused]] const auto lck = fs.lock();
  return unsafe_deserialize(offset, fs);
}

NormalizedExpectedValues NormalizedExpectedValues::unsafe_deserialize(std::streampos offset,
                                                                      filestream::FileStream &fs) {
  assert(offset >= 0);
  NormalizedExpectedValues nevs{};

  try {
    fs.unsafe_seekg(offset);
    const auto nNormExpectedValueVectors = static_cast<std::size_t>(fs.unsafe_read<std::int32_t>());
    for (std::size_t i = 0; i < nNormExpectedValueVectors; ++i) {
      nevs.emplace(NormalizedExpectedValuesBlock::unsafe_deserialize(fs.unsafe_tellg(), fs), true);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while deserializing a NormalizedExpectedValues object: {}"),
        e.what()));
  }

  return nevs;
}

NormalizationVectorIndexBlock::NormalizationVectorIndexBlock(
    std::string type_, std::uint32_t chrom_idx, std::string unit_, std::uint32_t bin_size,
    std::size_t position_, std::size_t n_bytes)
    : type(std::move(type_)),
      chrIdx(static_cast<std::int32_t>(chrom_idx)),
      unit(std::move(unit_)),
      binSize(static_cast<std::int32_t>(bin_size)),
      position(static_cast<std::int64_t>(position_)),
      nBytes(static_cast<std::int64_t>(n_bytes)) {}

bool NormalizationVectorIndexBlock::operator<(
    const NormalizationVectorIndexBlock &other) const noexcept {
  if (type != other.type) {
    return type < other.type;
  }
  if (chrIdx != other.chrIdx) {
    return chrIdx < other.chrIdx;
  }
  if (unit != other.unit) {
    return unit < other.unit;
  }
  return binSize < other.binSize;
}

std::string NormalizationVectorIndexBlock::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(type);
    buffer.write(chrIdx);
    buffer.write(unit);
    buffer.write(binSize);
    buffer.write(position);
    buffer.write(nBytes);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while serializing a NormalizationVectorIndexBlock object: {}"),
        e.what()));
  }

  return buffer.get();
}

NormalizationVectorIndexBlock NormalizationVectorIndexBlock::deserialize(
    std::streampos offset, filestream::FileStream &fs) {
  [[maybe_unused]] const auto lck = fs.lock();
  return unsafe_deserialize(offset, fs);
}

NormalizationVectorIndexBlock NormalizationVectorIndexBlock::unsafe_deserialize(
    std::streampos offset, filestream::FileStream &fs) {
  assert(offset >= 0);
  NormalizationVectorIndexBlock nvib{};

  try {
    fs.unsafe_seekg(offset);
    fs.unsafe_getline(nvib.type, '\0');
    nvib.chrIdx = fs.unsafe_read<std::int32_t>();
    fs.unsafe_getline(nvib.unit, '\0');
    nvib.binSize = fs.unsafe_read<std::int32_t>();
    nvib.position = fs.unsafe_read<std::int64_t>();
    nvib.nBytes = fs.unsafe_read<std::int64_t>();
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while deserializing a NormalizationVectorIndexBlock object: {}"),
        e.what()));
  }

  return nvib;
}

std::int32_t NormalizationVectorIndex::nNormVectors() const noexcept {
  return static_cast<std::int32_t>(_norm_vect_idx.size());
}

const std::vector<NormalizationVectorIndexBlock> &
NormalizationVectorIndex::normalizationVectorIndex() const noexcept {
  return _norm_vect_idx;
}

void NormalizationVectorIndex::emplace_back(NormalizationVectorIndexBlock blk) {
  _norm_vect_idx.emplace_back(std::move(blk));
}

std::string NormalizationVectorIndex::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  try {
    buffer.write(nNormVectors());

    for (const auto &nv : _norm_vect_idx) {
      std::ignore = nv.serialize(buffer, false);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a NormalizationVectorIndex object: {}"),
        e.what()));
  }

  return buffer.get();
}

NormalizationVectorIndex NormalizationVectorIndex::deserialize(std::streampos offset,
                                                               filestream::FileStream &fs) {
  [[maybe_unused]] const auto lck = fs.lock();
  return unsafe_deserialize(offset, fs);
}

NormalizationVectorIndex NormalizationVectorIndex::unsafe_deserialize(std::streampos offset,
                                                                      filestream::FileStream &fs) {
  assert(offset >= 0);
  NormalizationVectorIndex nvi{};

  try {
    fs.unsafe_seekg(offset);
    const auto nNormVectors = static_cast<std::size_t>(fs.unsafe_read<std::int32_t>());
    for (std::size_t i = 0; i < nNormVectors; ++i) {
      nvi.emplace_back(NormalizationVectorIndexBlock::unsafe_deserialize(fs.unsafe_tellg(), fs));
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while deserializing a NormalizationVectorIndex object: {}"),
        e.what()));
  }
  return nvi;
}

}  // namespace hictk::hic::internal
