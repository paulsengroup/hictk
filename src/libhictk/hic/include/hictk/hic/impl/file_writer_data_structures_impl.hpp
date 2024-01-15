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
inline bool MatrixInteractionBlock<N>::Pixel::operator<(const Pixel &other) const noexcept {
  return column < other.column;
}

template <typename N>
inline std::size_t MatrixInteractionBlock<N>::size() const noexcept {
  return static_cast<std::size_t>(nRecords);
}

template <typename N>
inline double MatrixInteractionBlock<N>::sum() const noexcept {
  return _sum;
}

template <typename N>
inline void MatrixInteractionBlock<N>::emplace_back(hictk::Pixel<N> &&p) {
  _sum += conditional_static_cast<double>(p.count);

  const auto row = static_cast<std::int32_t>(p.coords.bin2.rel_id());
  const auto col = static_cast<std::int32_t>(p.coords.bin1.rel_id());

  _min_col = std::min(col, _min_col);
  _max_col = std::max(col, _max_col);

  binRowOffset = std::min(binRowOffset, row);
  binColumnOffset = std::min(binColumnOffset, col);

  auto match1 = _interactions.find(row);
  if (match1 != _interactions.end()) {
    auto &pixels = match1->second;
    auto [it, inserted] = pixels.emplace(Pixel{col, p.count});
    nRecords += inserted;
    if (!inserted) {
      it->count += p.count;
    }
  } else {
    nRecords++;
    _interactions.emplace(row, Row{Pixel{col, p.count}});
  }
}

template <typename N>
inline void MatrixInteractionBlock<N>::finalize() {
  const auto size_lor = compute_size_lor_repr();
  const auto size_dense = compute_size_dense_repr();
  const auto width = compute_dense_width();

  const auto use_lor = size_lor < size_dense && width <= std::numeric_limits<std::int16_t>::max();

  useFloatContact = 1;
  useIntXPos = 1;
  useIntYPos = 1;
  matrixRepresentation = use_lor ? 1 : 2;

  // his can overflow, but it's ok because in this case use_lor=true
  w = static_cast<std::int16_t>(width);
}

template <typename N>
inline auto MatrixInteractionBlock<N>::operator()() const noexcept
    -> const phmap::btree_map<RowID, Row> & {
  return _interactions;
}

template <typename N>
inline std::string MatrixInteractionBlock<N>::serialize(BinaryBuffer &buffer,
                                                        libdeflate_compressor &compressor,
                                                        std::string &compression_buffer,
                                                        bool clear) const {
  if (matrixRepresentation == 1) {
    return serialize_lor(buffer, compressor, compression_buffer, clear);
  }
  return serialize_dense(buffer, compressor, compression_buffer, clear);
}

template <typename N>
inline std::size_t MatrixInteractionBlock<N>::compute_size_lor_repr() const noexcept {
  std::size_t size_ = sizeof(nRecords) + sizeof(binColumnOffset) + sizeof(binRowOffset) +
                      sizeof(useFloatContact) + sizeof(useIntXPos) + sizeof(useIntYPos) +
                      sizeof(matrixRepresentation);

  // compute space taken up by rows
  size_ += (_interactions.size() * sizeof(std::int32_t)) + sizeof(std::int32_t);

  // compute space taken up by columns
  size_ += size() * (sizeof(std::int32_t) + sizeof(N));

  return size_;
}

template <typename N>
inline std::size_t MatrixInteractionBlock<N>::compute_size_dense_repr() const noexcept {
  const auto width = compute_dense_width();
  const auto npixels = width * width;

  const std::size_t size_ = sizeof(nRecords) + sizeof(binColumnOffset) + sizeof(binRowOffset) +
                            sizeof(useFloatContact) + sizeof(useIntXPos) + sizeof(useIntYPos) +
                            sizeof(matrixRepresentation);
  return size_ + (sizeof(std::int32_t) + sizeof(std::int16_t)) + (npixels * sizeof(N));
}

template <typename N>
inline std::size_t MatrixInteractionBlock<N>::compute_dense_width() const noexcept {
  const auto min_row = _interactions.begin()->first;
  const auto max_row = (--_interactions.end())->first;
  const auto height = max_row - min_row;

  const auto width = _max_col - _min_col;

  return static_cast<std::size_t>(std::max(height, width) + 1);
}

template <typename N>
inline std::string MatrixInteractionBlock<N>::serialize_lor(BinaryBuffer &buffer,
                                                            libdeflate_compressor &compressor,
                                                            std::string &compression_buffer,
                                                            bool clear) const {
  assert(matrixRepresentation == 1);
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

  const auto rowCount = static_cast<std::int32_t>(_interactions.size());  // TODO support short
  buffer.write(rowCount);

  for (const auto &[row, pixels] : _interactions) {
    assert(static_cast<std::int32_t>(row) >= binRowOffset);
    const auto rowNumber = static_cast<std::int32_t>(row) - binRowOffset;  // TODO support short
    const auto recordCount = static_cast<std::int32_t>(pixels.size());     // TODO support short
    buffer.write(rowNumber);
    buffer.write(recordCount);

    assert(std::is_sorted(pixels.begin(), pixels.end()));
    for (const auto &[col, count] : pixels) {
      assert(col >= binColumnOffset);
      const auto binColumn = col - binColumnOffset;
      buffer.write(binColumn);
      buffer.write(count);
    }
  }

  compress(buffer.get(), compression_buffer, compressor);
  return compression_buffer;
}

template <typename N>
inline std::string MatrixInteractionBlock<N>::serialize_dense(BinaryBuffer &buffer,
                                                              libdeflate_compressor &compressor,
                                                              std::string &compression_buffer,
                                                              bool clear) const {
  assert(matrixRepresentation == 2);
  // TODO support representation using shorts

  if (clear) {
    buffer.clear();
  }

  const N fill_value = -32768;
  std::vector<N> counts(static_cast<std::size_t>(w) * static_cast<std::size_t>(w), fill_value);

  for (const auto &[row, pixels] : _interactions) {
    assert(row >= binRowOffset);
    const auto i = static_cast<std::size_t>(row - binRowOffset);
    for (const auto &[col, value] : pixels) {
      const auto j = static_cast<std::size_t>(col - binColumnOffset);
      const auto idx = (i * static_cast<std::size_t>(w)) + j;
      assert(idx < counts.size());
      counts[idx] = value;
    }
  }

  if constexpr (std::is_floating_point_v<N>) {
    std::transform(counts.begin(), counts.end(), counts.begin(), [&](const auto n) {
      return n == fill_value ? std::numeric_limits<N>::quiet_NaN() : n;
    });
  }

  buffer.write(nRecords);
  buffer.write(binColumnOffset);
  buffer.write(binRowOffset);
  buffer.write(useFloatContact);
  buffer.write(useIntXPos);
  buffer.write(useIntYPos);
  buffer.write(matrixRepresentation);

  buffer.write(static_cast<std::int32_t>(counts.size()));
  buffer.write(w);
  buffer.write(counts);

  compress(buffer.get(), compression_buffer, compressor);
  return compression_buffer;
}

template <typename N>
inline void MatrixInteractionBlock<N>::compress(const std::string &buffer_in,
                                                std::string &buffer_out,
                                                libdeflate_compressor &compressor) {
  assert(buffer_out.capacity() != 0);
  buffer_out.resize(buffer_out.capacity());
  while (true) {
    const auto compressed_size = libdeflate_zlib_compress(
        &compressor, buffer_in.data(), buffer_in.size(), buffer_out.data(), buffer_out.size());
    if (compressed_size != 0) {
      buffer_out.resize(compressed_size);
      break;
    }

    buffer_out.resize(buffer_out.size() * 2);
  }
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
