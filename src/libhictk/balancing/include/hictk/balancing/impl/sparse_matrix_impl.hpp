// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>
#include <zstd.h>

#include <cmath>
#include <iostream>
#include <iterator>
#include <nonstd/span.hpp>
#include <numeric>
#include <type_traits>

#include "hictk/cooler/cooler.hpp"
#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::balancing {

inline SparseMatrix::SparseMatrix(const hictk::BinTable& bins, std::uint32_t chrom_id)
    : _chrom_id(chrom_id == _gw_id ? 0 : chrom_id),
      _marg(chrom_id == _gw_id ? bins.size() : bins.subset(chrom_id).size()) {
  _chrom_offsets.push_back(0);
  for (const Chromosome& chrom : bins.chromosomes()) {
    const auto nbins = (chrom.size() + bins.bin_size() - 1) / bins.bin_size();
    _chrom_offsets.push_back(_chrom_offsets.back() + nbins);
  }
  _bin1_offsets.resize(_chrom_offsets.size(), 0);
}

inline bool SparseMatrix::empty() const noexcept { return size() == 0; }
inline std::size_t SparseMatrix::size() const noexcept { return _counts.size(); }

inline void SparseMatrix::clear() noexcept {
  _bin1_ids.clear();
  _bin2_ids.clear();
  _counts.clear();
}

inline void SparseMatrix::shrink_to_fit() noexcept {
  _bin1_ids.shrink_to_fit();
  _bin2_ids.shrink_to_fit();
  _counts.shrink_to_fit();
  _chrom_offsets.shrink_to_fit();
}

inline void SparseMatrix::finalize() {
  shrink_to_fit();
  for (std::size_t i = 1; i < _bin1_offsets.size(); ++i) {
    if (_bin1_offsets[i] == 0) {
      _bin1_offsets[i] = _bin1_offsets[i - 1];
    }
  }
}

inline const std::vector<std::size_t>& SparseMatrix::bin1_ids() const noexcept { return _bin1_ids; }
inline const std::vector<std::size_t>& SparseMatrix::bin2_ids() const noexcept { return _bin2_ids; }
inline const std::vector<double>& SparseMatrix::counts() const noexcept { return _counts; }
inline const std::vector<double>& SparseMatrix::margs() const noexcept { return _marg; }
inline const std::vector<std::size_t>& SparseMatrix::chrom_offsets() const noexcept {
  return _chrom_offsets;
}

inline void SparseMatrix::push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count) {
  if (!empty() && bin1_id >= _chrom_offsets[_chrom_id + 1]) {
    _chrom_id = static_cast<std::uint32_t>(
        std::lower_bound(_chrom_offsets.begin(), _chrom_offsets.end(), _bin1_ids.back()) -
        _chrom_offsets.begin());
    _bin1_offsets[_chrom_id] = size();
  }

  _bin1_ids.push_back(bin1_id);
  _bin2_ids.push_back(bin2_id);
  _counts.push_back(count);
}

inline SparseMatrixView SparseMatrix::subset(std::uint32_t chrom_id) const {
  assert(chrom_id + 1 < chrom_offsets().size());
  const auto i0 = _bin1_offsets[chrom_id];
  const auto i1 = _bin1_offsets[chrom_id + 1];

  const auto bin1_ids_ = nonstd::span(bin1_ids()).subspan(i0, i1 - i0);
  const auto bin2_ids_ = nonstd::span(bin2_ids()).subspan(i0, i1 - i0);
  const auto counts_ = nonstd::span(counts()).subspan(i0, i1 - i0);

  const auto j0 = chrom_offsets()[chrom_id];
  const auto j1 = chrom_offsets()[chrom_id + 1];

  return {bin1_ids_, bin2_ids_, counts_, j0, j1 - j0};
}

inline SparseMatrixView SparseMatrix::view() const {
  const auto bin1_ids_ = nonstd::span(bin1_ids());
  const auto bin2_ids_ = nonstd::span(bin2_ids());
  const auto counts_ = nonstd::span(counts());

  return {bin1_ids_, bin2_ids_, counts_, 0, _marg.size()};
}

void SparseMatrix::serialize(std::fstream& fs, int compression_lvl) const {
  const auto size_ = size();
  fs.write(reinterpret_cast<const char*>(&size_), sizeof(std::size_t));

  const auto tmpbuff_size = ZSTD_compressBound(size() * sizeof(std::uint64_t));
  std::string tmpbuff(tmpbuff_size, '\0');

  std::size_t compressed_size =
      ZSTD_compress(reinterpret_cast<void*>(tmpbuff.data()), tmpbuff.size() * sizeof(char),
                    reinterpret_cast<const void*>(_bin1_ids.data()), size() * sizeof(std::uint64_t),
                    compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(reinterpret_cast<const char*>(&compressed_size), sizeof(std::size_t));
  fs.write(tmpbuff.data(), static_cast<std::streamsize>(compressed_size));

  compressed_size =
      ZSTD_compress(reinterpret_cast<void*>(tmpbuff.data()), tmpbuff.size() * sizeof(char),
                    reinterpret_cast<const void*>(_bin2_ids.data()), size() * sizeof(std::uint64_t),
                    compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(reinterpret_cast<const char*>(&compressed_size), sizeof(std::size_t));
  fs.write(tmpbuff.data(), static_cast<std::streamsize>(compressed_size));

  compressed_size = ZSTD_compress(
      reinterpret_cast<void*>(tmpbuff.data()), tmpbuff.size() * sizeof(char),
      reinterpret_cast<const void*>(_counts.data()), size() * sizeof(double), compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(reinterpret_cast<const char*>(&compressed_size), sizeof(std::size_t));
  fs.write(tmpbuff.data(), static_cast<std::streamsize>(compressed_size));
}

void SparseMatrix::deserialize(std::fstream& fs) {
  std::size_t size{};
  fs.read(reinterpret_cast<char*>(&size), sizeof(std::size_t));

  _bin1_ids.resize(size);
  _bin2_ids.resize(size);
  _counts.resize(size);

  std::string tmpbuff{};
  std::size_t compressed_size{};
  fs.read(reinterpret_cast<char*>(&compressed_size), sizeof(std::size_t));

  tmpbuff.resize(compressed_size);
  fs.read(tmpbuff.data(), static_cast<std::streamsize>(tmpbuff.size() * sizeof(char)));
  std::size_t decompressed_size =
      ZSTD_decompress(reinterpret_cast<char*>(_bin1_ids.data()), size * sizeof(std::uint64_t),
                      tmpbuff.data(), tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }

  fs.read(reinterpret_cast<char*>(&compressed_size), sizeof(std::size_t));
  tmpbuff.resize(compressed_size);
  fs.read(tmpbuff.data(), static_cast<std::streamsize>(tmpbuff.size() * sizeof(char)));
  decompressed_size =
      ZSTD_decompress(reinterpret_cast<char*>(_bin2_ids.data()), size * sizeof(std::uint64_t),
                      tmpbuff.data(), tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }

  fs.read(reinterpret_cast<char*>(&compressed_size), sizeof(std::size_t));
  tmpbuff.resize(compressed_size);
  fs.read(tmpbuff.data(), static_cast<std::streamsize>(tmpbuff.size() * sizeof(char)));
  decompressed_size =
      ZSTD_decompress(reinterpret_cast<char*>(_counts.data()), size * sizeof(double),
                      tmpbuff.data(), tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }
}

inline SparseMatrixChunked::SparseMatrixChunked(const BinTable& bins,
                                                std::filesystem::path tmp_file,
                                                std::size_t chunk_size, int compression_lvl)
    : _matrix(bins),
      _path(std::move(tmp_file)),
      _marg(bins.size()),
      _chrom_offsets(_matrix.chrom_offsets()),
      _bin1_offsets(_chrom_offsets.size(), 0),
      _chunk_size(chunk_size),
      _compression_lvl(compression_lvl) {
  _fs.exceptions(std::ios::badbit);
  _fs.open(_path, std::ios::out);

  _chrom_index.emplace(0, std::make_pair(std::size_t{}, std::size_t{}));
}

inline SparseMatrixChunked::~SparseMatrixChunked() noexcept {
  try {
    if (!_path.empty() && std::filesystem::exists(_path)) {
      std::filesystem::remove(_path);
    }
  } catch (...) {
  }
}

inline bool SparseMatrixChunked::empty() const noexcept { return size() == 0; }
inline std::size_t SparseMatrixChunked::size() const noexcept { return _size; }

inline const std::vector<double>& SparseMatrixChunked::margs() const noexcept { return _marg; }
inline const std::vector<std::size_t>& SparseMatrixChunked::chrom_offsets() const noexcept {
  return _chrom_offsets;
}

inline void SparseMatrixChunked::push_back(std::uint64_t bin1_id, std::uint64_t bin2_id,
                                           double count) {
  const auto beginning_of_new_chromosome = bin1_id >= _chrom_offsets[_chrom_id + 1];

  if (beginning_of_new_chromosome) {
    finalize_chromosome(_chrom_id);
    _chrom_id = static_cast<std::uint32_t>(
        std::lower_bound(_chrom_offsets.begin(), _chrom_offsets.end(), _matrix.bin1_ids().back()) -
        _chrom_offsets.begin());
  }

  if (_matrix.size() == _chunk_size || beginning_of_new_chromosome) {
    write_chunk();
  }

  _matrix.push_back(bin1_id, bin2_id, count);
  ++_size;
}

inline void SparseMatrixChunked::finalize() {
  if (!_matrix.empty()) {
    write_chunk();
  }
  _fs.close();
  _fs.open(_path, std::ios::in);
}

inline void SparseMatrixChunked::finalize_chromosome(std::uint32_t chrom_id) {
  _bin1_offsets[chrom_id + 1] = size();
  auto [it, inserted] =
      _chrom_index.try_emplace(chrom_id, std::make_pair(_index.size(), _index.size()));
  if (!inserted) {
    it->second.second = _index.size() + 1;
  }
}

inline SparseMatrixChunkedView SparseMatrixChunked::view() const {
  return {_path, _index, 0, _marg.size()};
}

inline SparseMatrixChunkedView SparseMatrixChunked::subset(std::uint32_t chrom_id) const {
  auto it = _chrom_index.find(chrom_id);
  if (it == _chrom_index.end()) {
    return {};
  }
  const auto& [first_offset, last_offset] = it->second;
  const auto i0 = chrom_offsets()[chrom_id];
  const auto i1 = chrom_offsets()[chrom_id + 1];

  return {_path, nonstd::span(_index).subspan(first_offset, last_offset - first_offset), i0,
          i1 - i0};
}

inline void SparseMatrixChunked::read_chunk(std::size_t chunk_id, SparseMatrix& buffer) {
  assert(chunk_id < _index.size());
  const auto offset = _index[chunk_id];

  std::fstream fs;
  _fs.exceptions(std::ios::badbit);
  fs.open(_path, std::ios::in);
  fs.seekg(offset);

  buffer.deserialize(fs);
}

inline void SparseMatrixChunked::write_chunk() {
  assert(!_matrix.empty());
  _index.push_back(_fs.tellg());
  _matrix.serialize(_fs, _compression_lvl);
  _matrix.clear();
  _chrom_index.try_emplace(_chrom_id, std::make_pair(_index.size(), _index.size()));
}

inline SparseMatrixView::SparseMatrixView(nonstd::span<const std::size_t> bin1_ids_,
                                          nonstd::span<const std::size_t> bin2_ids_,
                                          nonstd::span<const double> counts_,
                                          std::size_t bin1_offset, std::size_t num_bins)
    : _marg(num_bins),
      _bin1_offset(bin1_offset),
      bin1_ids(bin1_ids_),
      bin2_ids(bin2_ids_),
      counts(counts_) {}

inline bool SparseMatrixView::empty() const noexcept { return size() == 0; }
inline std::size_t SparseMatrixView::size() const noexcept { return counts.size(); }

inline const std::vector<double>& SparseMatrixView::margs() const noexcept { return _marg; }

inline const std::vector<double>& SparseMatrixView::marginalize() const {
  std::fill(_marg.begin(), _marg.end(), 0);
  for (std::size_t i = 0; i < size(); ++i) {
    const auto i1 = bin1_ids[i] - _bin1_offset;
    const auto i2 = bin2_ids[i] - _bin1_offset;

    _marg[i1] += counts[i];
    _marg[i2] += counts[i];
  }

  return _marg;
}

inline const std::vector<double>& SparseMatrixView::marginalize_nnz() const {
  std::fill(_marg.begin(), _marg.end(), 0);

  for (std::size_t i = 0; i < counts.size(); ++i) {
    const auto i1 = bin1_ids[i] - _bin1_offset;
    const auto i2 = bin2_ids[i] - _bin1_offset;

    _marg[i1] += counts[i] != 0;
    _marg[i2] += counts[i] != 0;
  }

  return _marg;
}

inline const std::vector<double>& SparseMatrixView::times_outer_product_marg(
    nonstd::span<const double> biases, nonstd::span<const double> weights) const {
  assert(biases.size() == _marg.size());
  assert(biases.size() == weights.size() || weights.empty());

  std::fill(_marg.begin(), _marg.end(), 0);
  for (std::size_t i = 0; i < size(); ++i) {
    const auto i1 = bin1_ids[i] - _bin1_offset;
    const auto i2 = bin2_ids[i] - _bin1_offset;
    const auto w1 = weights.empty() ? 1 : weights[i1];
    const auto w2 = weights.empty() ? 1 : weights[i2];
    const auto count = counts[i] * (w1 * biases[i1]) * (w2 * biases[i2]);

    _marg[i1] += count;
    _marg[i2] += count;
  }
  return _marg;
}

inline SparseMatrixChunkedView::SparseMatrixChunkedView(const std::filesystem::path& path,
                                                        nonstd::span<const std::streamoff> index,
                                                        std::size_t bin1_offset,
                                                        std::size_t num_bins)
    : _fs(path, std::ios::in),
      _index(index.begin(), index.end()),
      _marg(num_bins),
      _bin1_offset(bin1_offset) {}

inline bool SparseMatrixChunkedView::empty() const noexcept { return _index.empty(); }

inline const std::vector<double>& SparseMatrixChunkedView::margs() const noexcept { return _marg; }

inline const std::vector<double>& SparseMatrixChunkedView::marginalize() const {
  std::fill(_marg.begin(), _marg.end(), 0);

  for (const auto offset : _index) {
    _fs.seekg(offset);
    _matrix.deserialize(_fs);

    for (std::size_t i = 0; i < _matrix.counts().size(); ++i) {
      const auto i1 = _matrix.bin1_ids()[i] - _bin1_offset;
      const auto i2 = _matrix.bin2_ids()[i] - _bin1_offset;

      _marg[i1] += _matrix.counts()[i];
      _marg[i2] += _matrix.counts()[i];
    }
    if (_fs.peek() && _fs.eof()) {
      break;
    }
  }

  return _marg;
}

inline const std::vector<double>& SparseMatrixChunkedView::marginalize_nnz() const {
  std::fill(_marg.begin(), _marg.end(), 0);

  for (const auto offset : _index) {
    _fs.seekg(offset);
    _matrix.deserialize(_fs);

    for (std::size_t i = 0; i < _matrix.counts().size(); ++i) {
      const auto i1 = _matrix.bin1_ids()[i] - _bin1_offset;
      const auto i2 = _matrix.bin2_ids()[i] - _bin1_offset;

      _marg[i1] += _matrix.counts()[i] != 0;
      _marg[i2] += _matrix.counts()[i] != 0;
    }
    if (_fs.peek() && _fs.eof()) {
      break;
    }
  }

  return _marg;
}

inline const std::vector<double>& SparseMatrixChunkedView::times_outer_product_marg(
    nonstd::span<const double> biases, nonstd::span<const double> weights) const {
  assert(biases.size() == _marg.size());
  assert(biases.size() == weights.size() || weights.empty());

  std::fill(_marg.begin(), _marg.end(), 0);

  for (const auto offset : _index) {
    _fs.seekg(offset);
    _matrix.deserialize(_fs);

    for (std::size_t i = 0; i < _matrix.counts().size(); ++i) {
      const auto i1 = _matrix.bin1_ids()[i] - _bin1_offset;
      const auto i2 = _matrix.bin2_ids()[i] - _bin1_offset;
      const auto w1 = weights.empty() ? 1 : weights[i1];
      const auto w2 = weights.empty() ? 1 : weights[i2];
      const auto count = _matrix.counts()[i] * (w1 * biases[i1]) * (w2 * biases[i2]);

      _marg[i1] += count;
      _marg[i2] += count;
    }
    if (_fs.peek() && _fs.eof()) {
      break;
    }
  }
  return _marg;
}

}  // namespace hictk::balancing
