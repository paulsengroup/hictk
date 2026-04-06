// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <libdeflate.h>
#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include "hictk/binary_buffer.hpp"
#include "hictk/fmt/pixel.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

template <typename It>
void MatrixResolutionMetadata::set_block_metadata(It first_block, It last_block) {
  _block_metadata.clear();
  std::copy(first_block, last_block, std::back_inserter(_block_metadata));
  blockCount = static_cast<std::int32_t>(_block_metadata.size());
}

template <typename N>
bool MatrixInteractionBlock<N>::Pixel::operator<(const Pixel &other) const noexcept {
  return column < other.column;
}

template <typename N>
std::size_t MatrixInteractionBlock<N>::size() const noexcept {
  return static_cast<std::size_t>(nRecords);
}

template <typename N>
double MatrixInteractionBlock<N>::sum() const noexcept {
  return _sum;
}

template <typename N>
void MatrixInteractionBlock<N>::emplace_back(const hictk::Pixel<N> &p,
                                             std::uint32_t bin_id_offset) {
  try {
    _sum += conditional_static_cast<double>(p.count);

    assert(p.coords.bin1.rel_id() >= bin_id_offset);
    assert(p.coords.bin2.rel_id() >= bin_id_offset);
    const auto row = static_cast<std::int32_t>(p.coords.bin2.rel_id() - bin_id_offset);
    const auto col = static_cast<std::int32_t>(p.coords.bin1.rel_id() - bin_id_offset);

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
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while adding pixel {} to a MatrixInteractionBlock object: {}"),
        p, e.what()));
  }
}

template <typename N>
void MatrixInteractionBlock<N>::finalize() {
  try {
    const auto size_lor = compute_size_lor_repr();
    const auto size_dense = compute_size_dense_repr();
    const auto width = compute_dense_width();

    const auto use_lor =
        (size_lor < size_dense) || (width > std::numeric_limits<std::int16_t>::max());

    useFloatContact = 1;
    useIntXPos = 1;
    useIntYPos = 1;
    matrixRepresentation = use_lor ? 1 : 2;

    // this can overflow, but it's ok because in this case use_lor=true
    w = static_cast<std::int16_t>(width);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while finalizing a MatrixInteractionBlock object: {}"),
        e.what()));
  }
}

template <typename N>
auto MatrixInteractionBlock<N>::operator()() const noexcept
    -> const phmap::btree_map<RowID, Row> & {
  return _interactions;
}

template <typename N>
std::string MatrixInteractionBlock<N>::serialize(BinaryBuffer &buffer,
                                                 libdeflate_compressor &compressor,
                                                 std::string &compression_buffer,
                                                 bool clear) const {
  if (matrixRepresentation == 1) {
    return serialize_lor(buffer, compressor, compression_buffer, clear);
  }
  return serialize_dense(buffer, compressor, compression_buffer, clear);
}

template <typename N>
std::size_t MatrixInteractionBlock<N>::compute_size_lor_repr() const noexcept {
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
std::size_t MatrixInteractionBlock<N>::compute_size_dense_repr() const noexcept {
  const auto width = compute_dense_width();
  const auto npixels = width * width;

  const std::size_t size_ = sizeof(nRecords) + sizeof(binColumnOffset) + sizeof(binRowOffset) +
                            sizeof(useFloatContact) + sizeof(useIntXPos) + sizeof(useIntYPos) +
                            sizeof(matrixRepresentation);
  return size_ + (sizeof(std::int32_t) + sizeof(std::int16_t)) + (npixels * sizeof(N));
}

template <typename N>
std::size_t MatrixInteractionBlock<N>::compute_dense_width() const noexcept {
  const auto min_row = _interactions.begin()->first;
  const auto max_row = (--_interactions.end())->first;
  const auto height = max_row - min_row;

  const auto width = _max_col - _min_col;

  assert(height >= 0);
  assert(width >= 0);
  return static_cast<std::size_t>(std::max(height, width)) + 1;
}

template <typename N>
std::string MatrixInteractionBlock<N>::serialize_lor(BinaryBuffer &buffer,
                                                     libdeflate_compressor &compressor,
                                                     std::string &compression_buffer,
                                                     bool clear) const {
  assert(matrixRepresentation == 1);
  // TODO support representation using shorts

  if (clear) {
    buffer.clear();
  }

  try {
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
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a MatrixInteractionBlock using the sparse "
                   "representation: {}"),
        e.what()));
  }

  try {
    compress(buffer.get(), compression_buffer, compressor);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while compressing a serialized object of "
                               "MatrixInteractionBlock type "
                               "(sparse representation): {}"),
                    e.what()));
  }
  return compression_buffer;
}

template <typename N>
std::string MatrixInteractionBlock<N>::serialize_dense(BinaryBuffer &buffer,
                                                       libdeflate_compressor &compressor,
                                                       std::string &compression_buffer,
                                                       bool clear) const {
  assert(matrixRepresentation == 2);
  // TODO support representation using shorts

  if (clear) {
    buffer.clear();
  }

  try {
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
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while serializing a MatrixInteractionBlock using the dense "
                   "representation: {}"),
        e.what()));
  }

  try {
    compress(buffer.get(), compression_buffer, compressor);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while compressing a serialized object of "
                               "MatrixInteractionBlock type "
                               "(dense representation): {}"),
                    e.what()));
  }

  return compression_buffer;
}

template <typename N>
void MatrixInteractionBlock<N>::compress(const std::string &buffer_in, std::string &buffer_out,
                                         libdeflate_compressor &compressor) {
  buffer_out.resize(buffer_out.capacity());
  if (buffer_out.empty()) {
    buffer_out.resize(libdeflate_deflate_compress_bound(&compressor, buffer_in.size()));
  }

  while (true) {
    const auto compressed_size = libdeflate_zlib_compress(
        &compressor, buffer_in.data(), buffer_in.size(), buffer_out.data(), buffer_out.size());
    if (compressed_size != 0) {
      buffer_out.resize(compressed_size);
      break;
    }

    const auto new_size = std::max(
        buffer_out.size() * 2, libdeflate_deflate_compress_bound(&compressor, buffer_in.size()));
    buffer_out.resize(new_size);
  }
}

}  // namespace hictk::hic::internal
