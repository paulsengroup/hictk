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

inline VectorOfAtomicDecimals::VectorOfAtomicDecimals(std::size_t size_, std::size_t decimals)
    : _margsi(size_), _margsd(size_), _cfx(static_cast<std::uint64_t>(std::pow(10, decimals - 1))) {
  fill(0);
}

inline VectorOfAtomicDecimals::VectorOfAtomicDecimals(const VectorOfAtomicDecimals& other)
    : _margsi(other.size()), _margsd(other.size()), _cfx(other._cfx) {
  for (std::size_t i = 0; i < size(); ++i) {
    _margsi[i] = other._margsi[i].load();
  }
}

inline VectorOfAtomicDecimals& VectorOfAtomicDecimals::operator=(
    const VectorOfAtomicDecimals& other) {
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

inline double VectorOfAtomicDecimals::operator[](std::size_t i) const noexcept {
  assert(i < size());
  return decode(_margsi[i].load());
}

inline void VectorOfAtomicDecimals::add(std::size_t i, double n) noexcept {
  assert(i < size());
  _margsi[i] += encode(n);
}

inline void VectorOfAtomicDecimals::set(std::size_t i, double n) noexcept {
  assert(i < size());
  _margsi[i] = encode(n);
}

inline void VectorOfAtomicDecimals::multiply(const std::vector<double>& v) noexcept {
  assert(size() == v.size());
  for (std::size_t i = 0; i < size(); ++i) {
    const auto n = decode(_margsi[i]);
    _margsi[i] = encode(n * v[i]);
  }
}

inline const std::vector<double>& VectorOfAtomicDecimals::operator()() const noexcept {
  assert(_margsi.size() == _margsd.size());
  for (std::size_t i = 0; i < size(); ++i) {
    _margsd[i] = (*this)[i];
  }
  return _margsd;
}

inline std::vector<double>& VectorOfAtomicDecimals::operator()() noexcept {
  assert(_margsi.size() == _margsd.size());
  for (std::size_t i = 0; i < size(); ++i) {
    _margsd[i] = (*this)[i];
  }
  return _margsd;
}

inline void VectorOfAtomicDecimals::fill(double value) noexcept {
  for (auto& n : _margsi) {
    n = encode(value);
  }
}

inline void VectorOfAtomicDecimals::resize(std::size_t size_) {
  if (size_ != size()) {
    _margsi = std::vector<N>(size_);
  }
}

inline std::size_t VectorOfAtomicDecimals::size() const noexcept { return _margsi.size(); }
inline bool VectorOfAtomicDecimals::empty() const noexcept { return size() == 0; }

inline auto VectorOfAtomicDecimals::encode(double n) const noexcept -> I {
  return static_cast<I>(n * static_cast<double>(_cfx));
}

inline double VectorOfAtomicDecimals::decode(I n) const noexcept {
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

inline void SparseMatrix::serialize(filestream::FileStream& fs, std::string& tmpbuff,
                                    ZSTD_CCtx& ctx, int compression_lvl) const {
  fs.write(size());

  const auto tmpbuff_size = ZSTD_compressBound(size() * sizeof(std::uint64_t));
  tmpbuff.resize(tmpbuff_size);

  std::size_t compressed_size = ZSTD_compressCCtx(&ctx, reinterpret_cast<void*>(tmpbuff.data()),
                                                  tmpbuff.size() * sizeof(char),
                                                  reinterpret_cast<const void*>(_bin1_ids.data()),
                                                  size() * sizeof(std::uint64_t), compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(compressed_size);
  fs.write(tmpbuff.data(), compressed_size);

  compressed_size = ZSTD_compressCCtx(&ctx, reinterpret_cast<void*>(tmpbuff.data()),
                                      tmpbuff.size() * sizeof(char),
                                      reinterpret_cast<const void*>(_bin2_ids.data()),
                                      size() * sizeof(std::uint64_t), compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(compressed_size);
  fs.write(tmpbuff.data(), compressed_size);

  compressed_size = ZSTD_compressCCtx(
      &ctx, reinterpret_cast<void*>(tmpbuff.data()), tmpbuff.size() * sizeof(char),
      reinterpret_cast<const void*>(_counts.data()), size() * sizeof(double), compression_lvl);
  if (ZSTD_isError(compressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(compressed_size));
  }

  fs.write(compressed_size);
  fs.write(tmpbuff.data(), compressed_size);

  fs.flush();
}

inline void SparseMatrix::deserialize(filestream::FileStream& fs, std::string& tmpbuff,
                                      ZSTD_DCtx& ctx) {
  const auto size_ = fs.read<std::size_t>();

  _bin1_ids.resize(size_);
  _bin2_ids.resize(size_);
  _counts.resize(size_);

  auto compressed_size = fs.read<std::size_t>();

  fs.read(tmpbuff, compressed_size);
  std::size_t decompressed_size = ZSTD_decompressDCtx(
      &ctx, reinterpret_cast<char*>(_bin1_ids.data()), _bin1_ids.size() * sizeof(std::uint64_t),
      tmpbuff.data(), tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }

  fs.read(compressed_size);
  fs.read(tmpbuff, compressed_size);
  decompressed_size = ZSTD_decompressDCtx(&ctx, reinterpret_cast<char*>(_bin2_ids.data()),
                                          _bin2_ids.size() * sizeof(std::uint64_t), tmpbuff.data(),
                                          tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }

  fs.read(compressed_size);
  fs.read(tmpbuff, compressed_size);
  decompressed_size = ZSTD_decompressDCtx(&ctx, reinterpret_cast<char*>(_counts.data()),
                                          _counts.size() * sizeof(double), tmpbuff.data(),
                                          tmpbuff.size() * sizeof(char));
  if (ZSTD_isError(decompressed_size)) {
    throw std::runtime_error(ZSTD_getErrorName(decompressed_size));
  }
}

inline void SparseMatrix::marginalize(VectorOfAtomicDecimals& marg, BS::thread_pool* tpool,
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

  tpool->detach_blocks(std::size_t(0), size(), marginalize_impl);
  tpool->wait();
}

inline void SparseMatrix::marginalize_nnz(VectorOfAtomicDecimals& marg, BS::thread_pool* tpool,
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

  tpool->detach_blocks(std::size_t(0), size(), marginalize_nnz_impl);
  tpool->wait();
}

inline void SparseMatrix::times_outer_product_marg(VectorOfAtomicDecimals& marg,
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

  tpool->detach_blocks(std::size_t(0), size(), times_outer_product_marg_impl);
  tpool->wait();
}

inline void SparseMatrix::multiply(VectorOfAtomicDecimals& buffer, nonstd::span<const double> cfx,
                                   BS::thread_pool* tpool, bool init_buffer) const {
  buffer.resize(cfx.size());

  if (init_buffer) {
    buffer.fill(0);
  }

  auto matrix_mult_impl = [&](std::size_t istart, std::size_t iend) {
    for (auto i = istart; i < iend; ++i) {
      const auto bin1_id = _bin1_ids[i];
      const auto bin2_id = _bin2_ids[i];
      const auto count = _counts[i];
      const auto w1 = cfx[bin1_id];
      const auto w2 = cfx[bin2_id];

      const auto f = bin1_id == bin2_id ? 0.5 : 1.0;
      buffer.add(bin1_id, count * f * w2);
      buffer.add(bin2_id, count * f * w1);
    }
  };

  if (size() < 1'000'000 || !tpool) {
    matrix_mult_impl(0, size());
    return;
  }

  tpool->detach_blocks(std::size_t(0), size(), matrix_mult_impl);
  tpool->wait();
}

inline double SparseMatrix::compute_scaling_factor_for_scale(
    const std::vector<double>& weights) const {
  if (empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double sum = 0.0;
  double norm_sum = 0.0;

  for (std::size_t i = 0; i < size(); ++i) {
    const auto bin1_id = bin1_ids()[i];
    const auto bin2_id = bin2_ids()[i];
    const auto count = counts()[i];

    const auto w1 = weights[bin1_id];
    const auto w2 = weights[bin2_id];

    if (!std::isnan(w1) && !std::isnan(w2)) {
      const auto cfx = bin1_id != bin2_id ? 2.0 : 1.0;
      sum += count * cfx;
      norm_sum += (count * cfx) / (w1 * w2);
    }
  }

  return std::sqrt(norm_sum / sum);
}

inline SparseMatrixChunked::SparseMatrixChunked(std::filesystem::path tmp_file,
                                                std::size_t chunk_size, int compression_lvl)
    : _path(std::move(tmp_file)),
      _fs(filestream::FileStream::create(_path.string())),
      _chunk_size(chunk_size),
      _compression_lvl(compression_lvl),
      _zstd_cctx(ZSTD_createCCtx()),
      _zstd_dctx(ZSTD_createDCtx()) {}

inline SparseMatrixChunked::~SparseMatrixChunked() noexcept {
  try {
    if (!_path.empty() && std::filesystem::exists(_path)) {
      _fs = filestream::FileStream{};
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
  _fs = filestream::FileStream{};
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
  _fs = filestream::FileStream(_path.string());
}

inline void SparseMatrixChunked::marginalize(VectorOfAtomicDecimals& marg, BS::thread_pool* tpool,
                                             bool init_buffer) const {
  auto marginalize_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    filestream::FileStream fs(_path.string());
    auto matrix = _matrix;
    std::string buff{};
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(static_cast<std::streamoff>(offset));
      matrix.deserialize(fs, buff, *zstd_dctx);
      matrix.marginalize(marg, nullptr, false);
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

    tpool->detach_task([&, ii = i0, jj = i1]() { marginalize_impl(ii, jj); });
  }
  tpool->wait();
}

inline void SparseMatrixChunked::marginalize_nnz(VectorOfAtomicDecimals& marg,
                                                 BS::thread_pool* tpool, bool init_buffer) const {
  auto marginalize_nnz_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    filestream::FileStream fs(_path.string());
    auto matrix = _matrix;
    std::string buff{};
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(static_cast<std::streamoff>(offset));
      matrix.deserialize(fs, buff, *zstd_dctx);
      matrix.marginalize_nnz(marg, nullptr, false);
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

    tpool->detach_task([&, ii = i0, jj = i1]() { marginalize_nnz_impl(ii, jj); });
  }
  tpool->wait();
}

inline void SparseMatrixChunked::times_outer_product_marg(VectorOfAtomicDecimals& marg,
                                                          nonstd::span<const double> biases,
                                                          nonstd::span<const double> weights,
                                                          BS::thread_pool* tpool,
                                                          bool init_buffer) const {
  auto times_outer_product_marg_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    filestream::FileStream fs(_path.string());
    auto matrix = _matrix;
    std::string buff{};
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(static_cast<std::streamoff>(offset));
      matrix.deserialize(fs, buff, *zstd_dctx);
      matrix.times_outer_product_marg(marg, biases, weights, nullptr, false);
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

    tpool->detach_task([&, ii = i0, jj = i1]() { times_outer_product_marg_impl(ii, jj); });
  }
  tpool->wait();
}

inline void SparseMatrixChunked::multiply(VectorOfAtomicDecimals& buffer,
                                          nonstd::span<const double> cfx, BS::thread_pool* tpool,
                                          bool init_buffer) const {
  auto times_outer_product_marg_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    filestream::FileStream fs(_path.string());
    auto matrix = _matrix;
    std::string buff{};
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(static_cast<std::streamoff>(offset));
      matrix.deserialize(fs, buff, *zstd_dctx);
      matrix.multiply(buffer, cfx, nullptr, false);
    }
  };

  buffer.resize(cfx.size());
  if (init_buffer) {
    buffer.fill(0);
  }

  if (_index.size() == 1 || !tpool) {
    times_outer_product_marg_impl(0, _index.size());
    return;
  }

  const auto offsets = compute_chunk_offsets(_index.size(), tpool->get_thread_count());

  for (std::size_t i = 1; i < offsets.size(); ++i) {
    const auto i0 = offsets[i - 1];
    const auto i1 = offsets[i];

    tpool->detach_task([&, ii = i0, jj = i1]() { times_outer_product_marg_impl(ii, jj); });
  }
  tpool->wait();
}

inline double SparseMatrixChunked::compute_scaling_factor_for_scale(
    const std::vector<double>& weights) const {
  if (empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double sum = 0.0;
  double norm_sum = 0.0;

  std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
  filestream::FileStream fs(_path.string());
  std::string buff{};

  for (const auto& offset : _index) {
    fs.seekg(static_cast<std::streamoff>(offset));
    _matrix.deserialize(fs, buff, *zstd_dctx);
    for (std::size_t i = 0; i < _matrix.size(); ++i) {
      const auto bin1_id = _matrix.bin1_ids()[i];
      const auto bin2_id = _matrix.bin2_ids()[i];
      const auto count = _matrix.counts()[i];

      const auto w1 = weights[bin1_id];
      const auto w2 = weights[bin2_id];

      assert(std::isfinite(w1));
      assert(std::isfinite(w2));

      if (!std::isnan(w1) && !std::isnan(w2)) {
        const auto cfx = bin1_id != bin2_id ? 2.0 : 1.0;
        sum += count * cfx;
        norm_sum += (count * cfx) / (w1 * w2);
      }
    }
  }

  return std::sqrt(norm_sum / sum);
}

inline void SparseMatrixChunked::write_chunk() {
  assert(!_matrix.empty());
  _index.push_back(_fs.tellp());
  _matrix.finalize();
  _matrix.serialize(_fs, _buff, *_zstd_cctx, _compression_lvl);
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
