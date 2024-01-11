// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <libdeflate.h>
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

#include "hictk/hic/binary_buffer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

inline std::string MatrixMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(chr1Idx);
  buffer.write(chr2Idx);
  buffer.write(nResolutions);

  return buffer.get();
}

inline std::string MatrixBlockMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(blockNumber);
  buffer.write(blockPosition);
  buffer.write(blockSizeBytes);

  return buffer.get();
}

inline bool MatrixBlockMetadata::operator<(const MatrixBlockMetadata &other) const noexcept {
  return blockNumber < other.blockNumber;
}

inline bool MatrixResolutionMetadata::operator<(
    const MatrixResolutionMetadata &other) const noexcept {
  if (unit != other.unit) {
    return unit < other.unit;
  }
  return binSize < other.binSize;
}

inline std::string MatrixResolutionMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  assert(!_block_metadata.empty());

  if (clear) {
    buffer.clear();
  }

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

  for (const auto &blk : _block_metadata) {
    std::ignore = blk.serialize(buffer, false);
  }

  return buffer.get();
}

template <typename It>
inline void MatrixResolutionMetadata::set_block_metadata(It first_block, It last_block) {
  _block_metadata.clear();
  std::copy(first_block, last_block, std::back_inserter(_block_metadata));
  blockCount = static_cast<std::int32_t>(_block_metadata.size());
}

inline std::string MatrixBodyMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  std::ignore = matrixMetadata.serialize(buffer, clear);
  for (const auto &metadata : resolutionMetadata) {
    std::ignore = metadata.serialize(buffer, false);
  }

  return buffer.get();
}

template <typename N>
inline void MatrixInteractionBlock<N>::emplace_back(Pixel<N> &&p) {
  nRecords++;

  const auto row = static_cast<std::int32_t>(p.coords.bin1.rel_id());
  const auto col = static_cast<std::int32_t>(p.coords.bin2.rel_id());

  binRowOffset = std::min(binRowOffset, row);
  binColumnOffset = std::min(binColumnOffset, col);

  auto it = _interactions.find(col);
  if (it != _interactions.end()) {
    it->second.push_back(std::move(p));
  } else {
    _interactions.emplace(col, std::vector<Pixel<float>>{std::move(p)});
  }
}

template <typename N>
inline void MatrixInteractionBlock<N>::finalize() {
  for (auto &[_, v] : _interactions) {
    std::sort(v.begin(), v.end());
  }

  // TODO tweak
  useFloatContact = 1;
  useIntXPos = 1;
  useIntYPos = 1;
  matrixRepresentation = 1;
}

template <typename N>
inline auto MatrixInteractionBlock<N>::operator()() const noexcept
    -> const phmap::btree_map<ColID, Col> & {
  return _interactions;
}

template <typename N>
inline std::string MatrixInteractionBlock<N>::serialize(BinaryBuffer &buffer,
                                                        libdeflate_compressor &compressor,
                                                        std::string &compression_buffer,
                                                        bool clear) const {
  // TODO support dense layout
  // TODO support representation using shorts

  if (clear) {
    buffer.clear();
  }

  buffer.write(nRecords);
  buffer.write(binColumnOffset);
  buffer.write(binRowOffset);
  buffer.write(useFloatContact);
  buffer.write(useIntXPos);
  buffer.write(useIntYPos);
  buffer.write(matrixRepresentation);

  const auto colCount = static_cast<std::int32_t>(_interactions.size());  // TODO support short
  buffer.write(colCount);

  for (const auto &[col, pixels] : _interactions) {
    assert(static_cast<std::int32_t>(col) >= binColumnOffset);
    const auto rowNumber = static_cast<std::int32_t>(col) - binColumnOffset;  // TODO support short
    const auto recordCount = static_cast<std::int32_t>(pixels.size());        // TODO support short
    buffer.write(rowNumber);
    buffer.write(recordCount);

    assert(std::is_sorted(pixels.begin(), pixels.end()));
    for (const auto &p : pixels) {
      const auto bin_id = static_cast<std::int32_t>(p.coords.bin1.rel_id());
      assert(bin_id >= binRowOffset);
      const auto binColumn = bin_id - binRowOffset;
      const auto value = p.count;
      buffer.write(binColumn);
      buffer.write(value);
    }
  }

  assert(compression_buffer.capacity() != 0);
  compression_buffer.resize(compression_buffer.capacity());
  while (true) {
    const auto compressed_size =
        libdeflate_zlib_compress(&compressor, buffer.get().data(), buffer.get().size(),
                                 compression_buffer.data(), compression_buffer.size());
    if (compressed_size != 0) {
      compression_buffer.resize(compressed_size);
      break;
    }

    compression_buffer.resize(compression_buffer.size() * 2);
  }

  return compression_buffer;
}

inline std::string MasterIndex::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(key);
  buffer.write(position);
  buffer.write(size);

  return buffer.get();
}

inline ExpectedValuesBlock::ExpectedValuesBlock(std::string_view unit_, std::uint32_t bin_size,
                                                const std::vector<double> &weights,
                                                const std::vector<std::uint32_t> &chrom_ids,
                                                const std::vector<double> &scale_factors)
    : unit(std::string{unit_}),
      binSize(static_cast<std::int32_t>(bin_size)),
      nValues(static_cast<std::int32_t>(weights.size())),
      value(weights.size()),
      nChrScaleFactors(static_cast<std::int32_t>(chrom_ids.size())),
      chrIndex(chrom_ids.size()),
      chrScaleFactor(chrom_ids.size()) {
  std::transform(weights.begin(), weights.end(), value.begin(),
                 [](const auto n) { return static_cast<float>(n); });
  std::transform(chrom_ids.begin(), chrom_ids.end(), chrIndex.begin(),
                 [](const auto n) { return static_cast<std::int32_t>(n); });
  std::transform(scale_factors.begin(), scale_factors.end(), chrScaleFactor.begin(),
                 [](const auto n) { return static_cast<float>(n); });
}

inline std::string ExpectedValuesBlock::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(unit);
  buffer.write(binSize);
  buffer.write(nValues);
  buffer.write(value);
  buffer.write(nChrScaleFactors);

  assert(chrIndex.size() == chrScaleFactor.size());
  for (std::size_t i = 0; i < chrIndex.size(); ++i) {
    buffer.write(chrIndex[i]);
    buffer.write(chrScaleFactor[i]);
  }

  return buffer.get();
}

inline std::string ExpectedValues::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(nExpectedValueVectors);

  if (nExpectedValueVectors == 0) {
    return buffer.get();
  }

  for (const auto &ev : expectedValues) {
    std::ignore = ev.serialize(buffer, false);
  }

  return buffer.get();
}

inline std::string NormalizedExpectedValues::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(nNormExpectedValueVectors);

  return buffer.get();
}

inline std::string NormalizationVectorIndex::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(nNormVectors);

  return buffer.get();
}

inline std::string NormalizationVectorArray::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(nValues);

  return buffer.get();
}

inline std::string FooterV5::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  return masterIndex.serialize(buffer);
}
}  // namespace hictk::hic::internal
