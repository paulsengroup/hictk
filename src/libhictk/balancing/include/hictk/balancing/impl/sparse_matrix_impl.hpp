// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <zstd.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <mutex>
#include <nonstd/span.hpp>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

namespace hictk::balancing {

inline MargsVector::MargsVector(std::size_t size_, std::size_t decimals)
    : _margsi(size_), _margsd(size_), _cfx(static_cast<std::uint64_t>(std::pow(10, decimals - 1))) {
  fill(0);
}

inline MargsVector::MargsVector(const MargsVector& other)
    : _margsi(other.size()), _margsd(other.size()), _cfx(other._cfx) {
  for (std::size_t i = 0; i < size(); ++i) {
    _margsi[i] = other._margsi[i].load();
  }
}

inline MargsVector& MargsVector::operator=(const MargsVector& other) {
  if (this == &other) {
    return *this;
  }
  _margsi = std::vector<N>(other.size());
  for (std::size_t i = 0; i < size(); ++i) {
    _margsi[i] = other._margsi[i].load();
  }
  _margsd = other._margsd;
  _cfx = other._cfx;

  return *this;
}

inline double MargsVector::operator[](std::size_t i) const noexcept {
  assert(i < size());
  return decode(_margsi[i].load());
}

inline void MargsVector::add(std::size_t i, double n) noexcept {
  assert(i < size());
  _margsi[i] += encode(n);
}

inline const std::vector<double>& MargsVector::operator()() const noexcept {
  assert(_margsi.size() == _margsd.size());
  for (std::size_t i = 0; i < size(); ++i) {
    _margsd[i] = (*this)[i];
  }
  return _margsd;
}

inline std::vector<double>& MargsVector::operator()() noexcept {
  assert(_margsi.size() == _margsd.size());
  for (std::size_t i = 0; i < size(); ++i) {
    _margsd[i] = (*this)[i];
  }
  return _margsd;
}

inline void MargsVector::fill(double value) noexcept {
  for (auto& n : _margsi) {
    n = encode(value);
  }
}

inline void MargsVector::resize(std::size_t size_) {
  if (size_ != size()) {
    _margsi = std::vector<N>(size_);
  }
}

inline std::size_t MargsVector::size() const noexcept { return _margsi.size(); }
inline bool MargsVector::empty() const noexcept { return size() == 0; }

inline auto MargsVector::encode(double n) const noexcept -> I {
  return static_cast<I>(n * static_cast<double>(_cfx));
}

inline double MargsVector::decode(I n) const noexcept {
  return static_cast<double>(n) / static_cast<double>(_cfx);
}

inline bool SparseMatrix::empty() const noexcept { return size() == 0; }
inline std::size_t SparseMatrix::size() const noexcept { return _counts.size(); }

inline void SparseMatrix::clear(bool shrink_to_fit_) noexcept {
  _bin1_ids.clear();
  _bin2_ids.clear();
  _counts.clear();
  if (shrink_to_fit_) {
    shrink_to_fit();
  }
}

inline void SparseMatrix::shrink_to_fit() noexcept {
  _bin1_ids.shrink_to_fit();
  _bin2_ids.shrink_to_fit();
  _counts.shrink_to_fit();
}

inline void SparseMatrix::finalize() { shrink_to_fit(); }

inline const std::vector<std::uint64_t>& SparseMatrix::bin1_ids() const noexcept {
  return _bin1_ids;
}
inline const std::vector<std::uint64_t>& SparseMatrix::bin2_ids() const noexcept {
  return _bin2_ids;
}
inline const std::vector<double>& SparseMatrix::counts() const noexcept { return _counts; }

inline void SparseMatrix::push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count,
                                    std::size_t bin_offset) {
  assert(bin1_id >= bin_offset);
  assert(bin2_id >= bin1_id);

  _bin1_ids.push_back(bin1_id - bin_offset);
  _bin2_ids.push_back(bin2_id - bin_offset);
  _counts.push_back(count);
}

inline void SparseMatrix::serialize(std::fstream& fs, ZSTD_CCtx& ctx, int compression_lvl) const {
  auto size_ = size();
  fs.write(reinterpret_cast<const char*>(&size_), sizeof(std::size_t));

  const auto tmpbuff_size = ZSTD_compressBound(size() * sizeof(std::uint64_t));
  std::string tmpbuff(tmpbuff_size, '\0');

  std::size_t compressed_size = ZSTD_compressCCtx(&ctx, reinterpret_cast<void*>(tmpbuff.data()),
                                                  tmpbuff.size() * sizeof(char),
                                                  reinterpret_cast<const void*>(_bin1_ids.data()),
                                                  size() * sizeof(std::uint64_t), compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(reinterpret_cast<const char*>(&compressed_size), sizeof(std::size_t));
  fs.write(tmpbuff.data(), static_cast<std::streamsize>(compressed_size));

  compressed_size = ZSTD_compressCCtx(&ctx, reinterpret_cast<void*>(tmpbuff.data()),
                                      tmpbuff.size() * sizeof(char),
                                      reinterpret_cast<const void*>(_bin2_ids.data()),
                                      size() * sizeof(std::uint64_t), compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(reinterpret_cast<const char*>(&compressed_size), sizeof(std::size_t));
  fs.write(tmpbuff.data(), static_cast<std::streamsize>(compressed_size));

  compressed_size = ZSTD_compressCCtx(
      &ctx, reinterpret_cast<void*>(tmpbuff.data()), tmpbuff.size() * sizeof(char),
      reinterpret_cast<const void*>(_counts.data()), size() * sizeof(double), compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(reinterpret_cast<const char*>(&compressed_size), sizeof(std::size_t));
  fs.write(tmpbuff.data(), static_cast<std::streamsize>(compressed_size));

  fs.flush();
}

inline void SparseMatrix::deserialize(std::fstream& fs, ZSTD_DCtx& ctx) {
  std::size_t size_{};
  fs.read(reinterpret_cast<char*>(&size_), sizeof(std::size_t));

  _bin1_ids.resize(size_);
  _bin2_ids.resize(size_);
  _counts.resize(size_);

  std::string tmpbuff{};
  std::size_t compressed_size{};
  fs.read(reinterpret_cast<char*>(&compressed_size), sizeof(std::size_t));

  tmpbuff.resize(compressed_size);
  fs.read(tmpbuff.data(), static_cast<std::streamsize>(tmpbuff.size() * sizeof(char)));
  std::size_t decompressed_size = ZSTD_decompressDCtx(
      &ctx, reinterpret_cast<char*>(_bin1_ids.data()), _bin1_ids.size() * sizeof(std::uint64_t),
      tmpbuff.data(), tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }

  fs.read(reinterpret_cast<char*>(&compressed_size), sizeof(std::size_t));
  tmpbuff.resize(compressed_size);
  fs.read(tmpbuff.data(), static_cast<std::streamsize>(tmpbuff.size() * sizeof(char)));
  decompressed_size = ZSTD_decompressDCtx(&ctx, reinterpret_cast<char*>(_bin2_ids.data()),
                                          _bin2_ids.size() * sizeof(std::uint64_t), tmpbuff.data(),
                                          tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }

  fs.read(reinterpret_cast<char*>(&compressed_size), sizeof(std::size_t));
  tmpbuff.resize(compressed_size);
  fs.read(tmpbuff.data(), static_cast<std::streamsize>(tmpbuff.size() * sizeof(char)));
  decompressed_size = ZSTD_decompressDCtx(&ctx, reinterpret_cast<char*>(_counts.data()),
                                          _counts.size() * sizeof(double), tmpbuff.data(),
                                          tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }
}

inline void SparseMatrix::marginalize(MargsVector& marg, BS::thread_pool* tpool,
                                      bool init_buffer) const {
  assert(!marg.empty());
  if (init_buffer) {
    marg.fill(0);
  }

  auto marginalize_impl = [&](std::size_t istart, std::size_t iend) {
    for (auto i = istart; i < iend; ++i) {
      const auto i1 = _bin1_ids[i];
      const auto i2 = _bin2_ids[i];

      if (_counts[i] != 0) {
        marg.add(i1, _counts[i]);
        marg.add(i2, _counts[i]);
      }
    }
  };

  if (size() < 1'000'000 || !tpool) {
    marginalize_impl(0, size());
    return;
  }

  tpool->push_loop(0, size(), marginalize_impl);
  tpool->wait_for_tasks();
}

inline void SparseMatrix::marginalize_nnz(MargsVector& marg, BS::thread_pool* tpool,
                                          bool init_buffer) const {
  if (init_buffer) {
    marg.fill(0);
  }

  auto marginalize_nnz_impl = [&](std::size_t istart, std::size_t iend) {
    for (auto i = istart; i < iend; ++i) {
      const auto i1 = _bin1_ids[i];
      const auto i2 = _bin2_ids[i];

      if (_counts[i] != 0) {
        marg.add(i1, _counts[i] != 0);
        marg.add(i2, _counts[i] != 0);
      }
    }
  };

  if (size() < 1'000'000 || !tpool) {
    marginalize_nnz_impl(0, size());
    return;
  }

  tpool->push_loop(0, size(), marginalize_nnz_impl);
  tpool->wait_for_tasks();
}

inline void SparseMatrix::times_outer_product_marg(MargsVector& marg,
                                                   nonstd::span<const double> biases,
                                                   nonstd::span<const double> weights,
                                                   BS::thread_pool* tpool, bool init_buffer) const {
  assert(biases.size() == weights.size() || weights.empty());
  marg.resize(biases.size());

  if (init_buffer) {
    marg.fill(0);
  }

  auto times_outer_product_marg_impl = [&](std::size_t istart, std::size_t iend) {
    for (auto i = istart; i < iend; ++i) {
      const auto i1 = _bin1_ids[i];
      const auto i2 = _bin2_ids[i];
      const auto w1 = weights.empty() ? 1 : weights[i1];
      const auto w2 = weights.empty() ? 1 : weights[i2];
      const auto count = _counts[i] * (w1 * biases[i1]) * (w2 * biases[i2]);

      if (count != 0) {
        marg.add(i1, count);
        marg.add(i2, count);
      }
    }
  };

  if (size() < 1'000'000 || !tpool) {
    times_outer_product_marg_impl(0, size());
    return;
  }

  tpool->push_loop(0, size(), times_outer_product_marg_impl);
  tpool->wait_for_tasks();
}

inline SparseMatrixChunked::SparseMatrixChunked(std::filesystem::path tmp_file,
                                                std::size_t chunk_size, int compression_lvl)
    : _path(std::move(tmp_file)),
      _chunk_size(chunk_size),
      _compression_lvl(compression_lvl),
      _zstd_cctx(ZSTD_createCCtx()),
      _zstd_dctx(ZSTD_createDCtx()) {
  _fs.exceptions(std::ios::badbit);
  _fs.open(_path, std::ios::out | std::ios::binary);
}

inline SparseMatrixChunked::~SparseMatrixChunked() noexcept {
  try {
    if (!_path.empty() && std::filesystem::exists(_path)) {
      std::filesystem::remove(_path);
    }
    // NOLINTNEXTLINE
  } catch (...) {
  }
}

inline bool SparseMatrixChunked::empty() const noexcept { return size() == 0; }
inline std::size_t SparseMatrixChunked::size() const noexcept { return _size; }
inline void SparseMatrixChunked::clear(bool shrink_to_fit_) {
  _index.clear();
  _fs.close();
  std::filesystem::remove(_path);
  _path = "";
  _size = 0;
  _matrix.clear(shrink_to_fit_);
}

inline void SparseMatrixChunked::push_back(std::uint64_t bin1_id, std::uint64_t bin2_id,
                                           double count, std::size_t bin_offset) {
  if (_matrix.size() == _chunk_size) {
    write_chunk();
  }

  _matrix.push_back(bin1_id, bin2_id, count, bin_offset);
  ++_size;
}

inline void SparseMatrixChunked::finalize() {
  if (!_matrix.empty()) {
    write_chunk();
  }
  _fs.open(_path, std::ios::in | std::ios::binary);
}

inline void SparseMatrixChunked::marginalize(MargsVector& marg, BS::thread_pool* tpool,
                                             bool init_buffer) const {
  auto marginalize_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    std::fstream fs{};
    fs.exceptions(_fs.exceptions());
    fs.open(_path, std::ios::in | std::ios::binary);
    auto matrix = _matrix;
    MargsVector marg_local(marg.size());
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(offset);
      matrix.deserialize(fs, *zstd_dctx);
      matrix.marginalize(marg_local, nullptr, false);
    }

    for (std::size_t i = 0; i < marg_local.size(); ++i) {
      if (marg_local[i] != 0) {
        marg.add(i, marg_local[i]);
      }
    }
  };

  assert(!marg.empty());
  if (init_buffer) {
    marg.fill(0);
  }

  if (_index.size() == 1 || !tpool) {
    marginalize_impl(0, _index.size());
    return;
  }

  const auto offsets = compute_chunk_offsets(_index.size(), tpool->get_thread_count());

  for (std::size_t i = 1; i < offsets.size(); ++i) {
    const auto i0 = offsets[i - 1];
    const auto i1 = offsets[i];

    tpool->push_task(marginalize_impl, i0, i1);
  }
  tpool->wait_for_tasks();
}

inline void SparseMatrixChunked::marginalize_nnz(MargsVector& marg, BS::thread_pool* tpool,
                                                 bool init_buffer) const {
  auto marginalize_nnz_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    std::fstream fs{};
    fs.exceptions(_fs.exceptions());
    fs.open(_path, std::ios::in | std::ios::binary);
    auto matrix = _matrix;
    MargsVector marg_local(marg.size());
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(offset);
      matrix.deserialize(fs, *zstd_dctx);
      matrix.marginalize_nnz(marg, nullptr, false);
    }
    for (std::size_t i = 0; i < marg_local.size(); ++i) {
      if (marg_local[i] != 0) {
        marg.add(i, marg_local[i]);
      }
    }
  };

  assert(!marg.empty());
  if (init_buffer) {
    marg.fill(0);
  }

  if (_index.size() == 1 || !tpool) {
    marginalize_nnz_impl(0, _index.size());
    return;
  }
  const auto offsets = compute_chunk_offsets(_index.size(), tpool->get_thread_count());

  for (std::size_t i = 1; i < offsets.size(); ++i) {
    const auto i0 = offsets[i - 1];
    const auto i1 = offsets[i];

    tpool->push_task(marginalize_nnz_impl, i0, i1);
  }
  tpool->wait_for_tasks();
}

inline void SparseMatrixChunked::times_outer_product_marg(MargsVector& marg,
                                                          nonstd::span<const double> biases,
                                                          nonstd::span<const double> weights,
                                                          BS::thread_pool* tpool,
                                                          bool init_buffer) const {
  auto times_outer_product_marg_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    std::fstream fs{};
    fs.exceptions(_fs.exceptions());
    fs.open(_path, std::ios::in | std::ios::binary);
    auto matrix = _matrix;
    MargsVector marg_local(marg.size());
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(offset);
      matrix.deserialize(fs, *zstd_dctx);
      matrix.times_outer_product_marg(marg_local, biases, weights, nullptr, false);
    }
    for (std::size_t i = 0; i < marg.size(); ++i) {
      if (marg_local[i] != 0) {
        marg.add(i, marg_local[i]);
      }
    }
  };

  assert(biases.size() == weights.size() || weights.empty());
  marg.resize(biases.size());
  if (init_buffer) {
    marg.fill(0);
  }

  if (_index.size() == 1 || !tpool) {
    times_outer_product_marg_impl(0, _index.size());
    return;
  }

  const auto offsets = compute_chunk_offsets(_index.size(), tpool->get_thread_count());

  for (std::size_t i = 1; i < offsets.size(); ++i) {
    const auto i0 = offsets[i - 1];
    const auto i1 = offsets[i];

    tpool->push_task(times_outer_product_marg_impl, i0, i1);
  }
  tpool->wait_for_tasks();
}

inline void SparseMatrixChunked::write_chunk() {
  assert(!_matrix.empty());
  _index.push_back(_fs.tellg());
  _matrix.finalize();
  _matrix.serialize(_fs, *_zstd_cctx, _compression_lvl);
  _matrix.clear();
}

inline std::vector<std::size_t> SparseMatrixChunked::compute_chunk_offsets(std::size_t size,
                                                                           std::size_t num_chunks) {
  std::vector<std::size_t> offsets{};
  if (size < num_chunks) {
    offsets.resize(size + 1, 1);
    offsets.front() = 0;
  } else {
    const auto n = size / num_chunks;
    offsets.resize(num_chunks + 1, n);
    offsets.front() = 0;
    auto tot = n * num_chunks;

    for (std::size_t i = 1; i < offsets.size(); ++i) {
      if (tot == size) {
        break;
      }
      offsets[i]++;
      tot++;
    }
  }
  for (std::size_t i = 1; i < offsets.size(); ++i) {
    offsets[i] += offsets[i - 1];
  }

  return offsets;
}

}  // namespace hictk::balancing
