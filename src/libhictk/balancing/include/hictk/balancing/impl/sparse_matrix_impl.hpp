// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <zstd.h>

#include <BS_thread_pool.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
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

#include "hictk/common.hpp"

namespace hictk::balancing::internal {

inline AtomicBitSet::AtomicBitSet(std::size_t size_, bool value)
    : _buff(size_ / sizeof(I)), _size(size_) {
  fill(value);
}

inline AtomicBitSet::AtomicBitSet(const AtomicBitSet& other)
    : _buff(other._buff.size()), _size(other._size) {
  std::transform(other._buff.begin(), other._buff.end(), _buff.begin(),
                 [](const auto& n) { return n.load(); });
}

inline AtomicBitSet& AtomicBitSet::operator=(const AtomicBitSet& other) {
  if (this == &other) {
    return *this;
  }
  const auto smallest_size = static_cast<std::ptrdiff_t>(std::min(size(), other.size()));
  _buff = std::vector<std::atomic<I>>(other._buff.size());
  _size = other._size;

  std::transform(other._buff.begin(), other._buff.begin() + smallest_size, _buff.begin(),
                 [](const auto& n) { return n.load(); });
  std::fill(_buff.begin() + smallest_size, _buff.end(), false);

  return *this;
}

inline void AtomicBitSet::atomic_set(std::size_t i, bool value) noexcept {
  assert(i < _size);
  const auto uint_offset = compute_offset(i);
  const auto bit_offset = i - uint_offset;

  if (HICTK_LIKELY(value)) {
    const auto byte = static_cast<std::uint8_t>(1U << bit_offset);
    _buff[uint_offset].fetch_or(byte);
  } else {
    const auto byte = static_cast<std::uint8_t>(~(1U << bit_offset));
    _buff[uint_offset].fetch_and(byte);
  }
}

inline bool AtomicBitSet::atomic_test(std::size_t i) const noexcept {
  assert(i < _size);

  return atomic_test(_buff, i);
}

inline std::size_t AtomicBitSet::size() const noexcept { return _size; }

inline bool AtomicBitSet::empty() const noexcept { return size() == 0; }

inline void AtomicBitSet::fill(bool value) noexcept {
  if (value) {
    std::fill(_buff.begin(), _buff.end(), std::numeric_limits<std::uint8_t>::max());
  } else {
    std::fill(_buff.begin(), _buff.end(), std::uint8_t{});
  }
}

inline void AtomicBitSet::resize(std::size_t size_, bool value) {
  if (size_ != size()) {
    const auto old_size = size();
    auto new_v = std::vector<std::atomic<I>>(size_ / sizeof(I));
    std::swap(new_v, _buff);
    _size = size_;

    for (std::size_t i = 0; i < std::min(old_size, size()); ++i) {
      atomic_set(i, atomic_test(new_v, i));
    }

    for (std::size_t i = std::min(old_size, size()); i < size(); ++i) {
      atomic_set(i, value);
    }
  }
}

inline std::size_t AtomicBitSet::compute_offset(std::size_t i) noexcept { return i / sizeof(I); }

inline bool AtomicBitSet::atomic_test(const std::vector<std::atomic<I>>& buff,
                                      std::size_t i) noexcept {
  assert(i / sizeof(I) < buff.size());

  const auto uint_offset = compute_offset(i);
  const auto bit_offset = i - uint_offset;

  const auto byte = buff[uint_offset].load();
  return byte & (1U << bit_offset);
}

inline VectorOfAtomicDecimals::VectorOfAtomicDecimals(std::size_t size_, std::uint64_t decimal_bits)
    : _margsi(size_),
      _margsd(size_),
      _nanmask(size_),
      _infmask(size_),
      _cfxi(2ULL << (decimal_bits - 1)),
      _cfxd(static_cast<double>(_cfxi)),
      _max_value(compute_max_value(static_cast<std::uint8_t>(decimal_bits))) {
  if (decimal_bits == 0 || decimal_bits > 63) {
    throw std::invalid_argument("decimal bits should be between 1 and 64");
  }

  fill(0);
}

inline VectorOfAtomicDecimals::VectorOfAtomicDecimals(const VectorOfAtomicDecimals& other)
    : _margsi(other.size()),
      _margsd(other.size()),
      _nanmask(other._nanmask),
      _infmask(other._infmask),
      _cfxi(other._cfxi),
      _cfxd(other._cfxd),
      _max_value(other._max_value) {
  std::transform(other._margsi.begin(), other._margsi.end(), _margsi.begin(),
                 [&](const auto& n) { return n.load(); });
}

inline VectorOfAtomicDecimals& VectorOfAtomicDecimals::operator=(
    const VectorOfAtomicDecimals& other) {
  if (this == &other) {
    return *this;
  }
  _margsi = std::vector<N>(other.size());
  std::transform(other._margsi.begin(), other._margsi.end(), _margsi.begin(),
                 [&](const auto& n) { return n.load(); });
  _margsd = other._margsd;

  _nanmask = other._nanmask;
  _infmask = other._infmask;
  _cfxi = other._cfxi;
  _cfxd = other._cfxd;
  _max_value = other._max_value;

  return *this;
}

inline double VectorOfAtomicDecimals::operator[](std::size_t i) const noexcept {
  assert(i < size());
  if (HICTK_UNLIKELY(_nanmask.atomic_test(i))) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (HICTK_UNLIKELY(_infmask.atomic_test(i))) {
    return std::numeric_limits<double>::infinity();
  }

  return decode(_margsi[i].load());
}

inline void VectorOfAtomicDecimals::atomic_add(std::size_t i, double n) noexcept {
  assert(i < size());

  if (HICTK_UNLIKELY(std::isnan(n))) {
    _nanmask.atomic_set(i, true);
    return;
  }

  if (HICTK_UNLIKELY(overflows(n))) {
    _infmask.atomic_set(i, true);
    return;
  }

  const auto old_value = _margsi[i].load();
  const auto en = encode(n);
  constexpr auto max_value = std::numeric_limits<std::uint64_t>::max();
  if (HICTK_UNLIKELY(max_value - en < old_value)) {
    _infmask.atomic_set(i, true);
    return;
  }

  const auto new_value = _margsi[i] += en;

  if (HICTK_UNLIKELY(new_value < old_value)) {
    _infmask.atomic_set(i, true);
  }
}

inline void VectorOfAtomicDecimals::set(std::size_t i, double n) noexcept {
  assert(i < size());

  if (HICTK_UNLIKELY(std::isnan(n))) {
    _nanmask.atomic_set(i, true);
    return;
  }

  if (HICTK_UNLIKELY(overflows(n))) {
    _infmask.atomic_set(i, true);
    _nanmask.atomic_set(i, false);
    return;
  }

  _margsi[i] = encode(n);
  _nanmask.atomic_set(i, false);
  _infmask.atomic_set(i, false);
}

inline void VectorOfAtomicDecimals::multiply(const std::vector<double>& v) noexcept {
  assert(size() == v.size());
  for (std::size_t i = 0; i < size(); ++i) {
    const auto n = decode(_margsi[i]) * v[i];

    if (HICTK_UNLIKELY(std::isnan(n))) {
      _nanmask.atomic_set(i, true);
      continue;
    }

    if (HICTK_UNLIKELY(overflows(n))) {
      _infmask.atomic_set(i, true);
      continue;
    }

    _margsi[i] = encode(n);
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

  _nanmask.fill(false);
  _infmask.fill(false);
}

inline void VectorOfAtomicDecimals::resize(std::size_t size_, double value) {
  if (size_ != size()) {
    auto new_v = std::vector<N>(size_);
    const auto i = static_cast<std::ptrdiff_t>(std::min(size(), size_));
    std::transform(_margsi.begin(), _margsi.begin() + i, new_v.begin(),
                   [&](const auto& n) { return n.load(); });

    const auto en = std::isfinite(value) && !overflows(value) ? encode(value) : std::uint64_t{};
    std::generate(_margsi.begin() + i, _margsi.end(), [&]() { return en; });
    std::swap(new_v, _margsi);
    _nanmask.resize(size_, std::isnan(value));
    _infmask.resize(size_, !std::isnan(value) && overflows(value));
  }
}

inline std::uint8_t VectorOfAtomicDecimals::decimal_bits() const noexcept {
#ifdef _MSC_VER
  return static_cast<std::uint8_t>(_tzcnt_u64(_cfxi));
#else
  return static_cast<std::uint8_t>(__builtin_ctzll(_cfxi));
#endif
}

inline std::size_t VectorOfAtomicDecimals::size() const noexcept { return _margsi.size(); }

inline bool VectorOfAtomicDecimals::empty() const noexcept { return size() == 0; }

inline std::pair<double, double> VectorOfAtomicDecimals::domain(bool include_inf) const noexcept {
  if (include_inf) {
    return std::make_pair(0.0, std::numeric_limits<double>::infinity());
  }
  return std::make_pair(0.0, _max_value);
}

inline auto VectorOfAtomicDecimals::encode(double n) const noexcept -> I {
  assert(std::isfinite(n));
  assert(n <= _max_value);

  const auto encoded_n = n * _cfxd;
  assert(encoded_n <= static_cast<double>(std::numeric_limits<std::uint64_t>::max()));
  return static_cast<I>(encoded_n);
}

inline double VectorOfAtomicDecimals::decode(I n) const noexcept {
  return static_cast<double>(n) / _cfxd;
}

constexpr bool VectorOfAtomicDecimals::overflows(double n) const noexcept { return n > _max_value; }

inline double VectorOfAtomicDecimals::compute_max_value(std::uint8_t decimal_bits) const noexcept {
  assert(decimal_bits < 64);
  return std::nextafter(static_cast<double>(std::numeric_limits<I>::max() >> decimal_bits), 0.0);
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

inline void SparseMatrix::reserve(std::size_t capacity) {
  _bin1_ids.reserve(capacity);
  _bin2_ids.reserve(capacity);
  _counts.reserve(capacity);
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

inline void SparseMatrix::marginalize(VectorOfAtomicDecimals& marg, bool init_buffer) const {
  assert(!marg.empty());
  if (init_buffer) {
    marg.fill(0);
  }

  for (std::size_t i = 0; i < size(); ++i) {
    const auto i1 = _bin1_ids[i];
    const auto i2 = _bin2_ids[i];

    if (_counts[i] != 0) {
      marg.atomic_add(i1, _counts[i]);
      marg.atomic_add(i2, _counts[i]);
    }
  }
}

inline void SparseMatrix::marginalize_nnz(VectorOfAtomicDecimals& marg, bool init_buffer) const {
  if (init_buffer) {
    marg.fill(0);
  }

  for (std::size_t i = 0; i < size(); ++i) {
    const auto i1 = _bin1_ids[i];
    const auto i2 = _bin2_ids[i];

    if (_counts[i] != 0) {
      marg.atomic_add(i1, _counts[i] != 0);
      marg.atomic_add(i2, _counts[i] != 0);
    }
  }
}

inline void SparseMatrix::times_outer_product_marg(VectorOfAtomicDecimals& marg,
                                                   nonstd::span<const double> biases,
                                                   nonstd::span<const double> weights,
                                                   bool init_buffer) const {
  assert(biases.size() == weights.size() || weights.empty());
  marg.resize(biases.size());

  if (init_buffer) {
    marg.fill(0);
  }

  for (std::size_t i = 0; i < size(); ++i) {
    const auto i1 = _bin1_ids[i];
    const auto i2 = _bin2_ids[i];
    const auto w1 = weights.empty() ? 1 : weights[i1];
    const auto w2 = weights.empty() ? 1 : weights[i2];
    const auto count = _counts[i] * (w1 * biases[i1]) * (w2 * biases[i2]);

    if (count != 0) {
      marg.atomic_add(i1, count);
      marg.atomic_add(i2, count);
    }
  }
}

inline void SparseMatrix::multiply(VectorOfAtomicDecimals& buffer, nonstd::span<const double> cfx,
                                   bool init_buffer) const {
  buffer.resize(cfx.size());

  if (init_buffer) {
    buffer.fill(0);
  }

  for (std::size_t i = 0; i < size(); ++i) {
    const auto bin1_id = _bin1_ids[i];
    const auto bin2_id = _bin2_ids[i];
    const auto count = _counts[i];
    const auto w1 = cfx[bin1_id];
    const auto w2 = cfx[bin2_id];

    const auto f = bin1_id == bin2_id ? 0.5 : 1.0;
    buffer.atomic_add(bin1_id, count * f * w2);
    buffer.atomic_add(bin2_id, count * f * w1);
  }
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

inline SparseMatrixChunked::SparseMatrixChunked(std::size_t chunk_size)
    : _chunks(1), _chunk_size(chunk_size) {
  assert(chunk_size != 0);
  _chunks.back().reserve(_chunk_size);
}

inline bool SparseMatrixChunked::empty() const noexcept { return size() == 0; }

inline std::size_t SparseMatrixChunked::size() const noexcept { return _size; }

inline std::size_t SparseMatrixChunked::num_chunks() const noexcept {
  return empty() ? std::size_t{} : _chunks.size();
}

inline std::size_t SparseMatrixChunked::chunk_size() const noexcept { return _chunk_size; }

inline void SparseMatrixChunked::shrink_to_fit() noexcept {
  assert(!_chunks.empty());
  _chunks.back().shrink_to_fit();
}

inline void SparseMatrixChunked::clear(bool shrink_to_fit_) {
  _chunks.resize(1);
  _chunks.back().clear();
  _size = 0;

  if (shrink_to_fit_) {
    shrink_to_fit();
  }
}

inline void SparseMatrixChunked::push_back(std::uint64_t bin1_id, std::uint64_t bin2_id,
                                           double count, std::size_t bin_offset) {
  assert(!_chunks.empty());
  if (HICTK_UNLIKELY(_chunks.back().size() == _chunk_size)) {
    auto& chunk = _chunks.emplace_back(SparseMatrix{});
    chunk.reserve(_chunk_size);
  }

  _chunks.back().push_back(bin1_id, bin2_id, count, bin_offset);
  ++_size;
}

inline void SparseMatrixChunked::finalize() { shrink_to_fit(); }

inline void SparseMatrixChunked::marginalize(VectorOfAtomicDecimals& marg, BS::thread_pool* tpool,
                                             bool init_buffer) const {
  auto marginalize_impl = [&](std::size_t istart, std::size_t iend) {
    for (std::size_t i = istart; i < iend; ++i) {
      _chunks[i].marginalize(marg, false);
    }
  };

  assert(!marg.empty());
  if (init_buffer) {
    marg.fill(0);
  }

  if (num_chunks() < 2 || !tpool) {
    marginalize_impl(0, num_chunks());
    return;
  }

  tpool->detach_blocks(std::size_t(0), num_chunks(), marginalize_impl);
  tpool->wait();
}

inline void SparseMatrixChunked::marginalize_nnz(VectorOfAtomicDecimals& marg,
                                                 BS::thread_pool* tpool, bool init_buffer) const {
  auto marginalize_nnz_impl = [&](std::size_t istart, std::size_t iend) {
    for (std::size_t i = istart; i < iend; ++i) {
      _chunks[i].marginalize_nnz(marg, false);
    }
  };

  assert(!marg.empty());
  if (init_buffer) {
    marg.fill(0);
  }

  if (num_chunks() < 2 || !tpool) {
    marginalize_nnz_impl(0, num_chunks());
    return;
  }

  tpool->detach_blocks(std::size_t(0), num_chunks(), marginalize_nnz_impl);
  tpool->wait();
}

inline void SparseMatrixChunked::times_outer_product_marg(VectorOfAtomicDecimals& marg,
                                                          nonstd::span<const double> biases,
                                                          nonstd::span<const double> weights,
                                                          BS::thread_pool* tpool,
                                                          bool init_buffer) const {
  auto times_outer_product_marg_impl = [&](std::size_t istart, std::size_t iend) {
    for (std::size_t i = istart; i < iend; ++i) {
      _chunks[i].times_outer_product_marg(marg, biases, weights, false);
    }
  };

  assert(biases.size() == weights.size() || weights.empty());
  marg.resize(biases.size());
  if (init_buffer) {
    marg.fill(0);
  }

  if (num_chunks() < 2 || !tpool) {
    times_outer_product_marg_impl(0, num_chunks());
    return;
  }

  tpool->detach_blocks(std::size_t(0), num_chunks(), times_outer_product_marg_impl);
  tpool->wait();
}

inline void SparseMatrixChunked::multiply(VectorOfAtomicDecimals& buffer,
                                          nonstd::span<const double> cfx, BS::thread_pool* tpool,
                                          bool init_buffer) const {
  auto multiply_impl = [&](std::size_t istart, std::size_t iend) {
    for (std::size_t i = istart; i < iend; ++i) {
      _chunks[i].multiply(buffer, cfx, false);
    }
  };

  buffer.resize(cfx.size());
  if (init_buffer) {
    buffer.fill(0);
  }

  if (num_chunks() < 2 || !tpool) {
    multiply_impl(0, num_chunks());
    return;
  }

  tpool->detach_blocks(std::size_t(0), num_chunks(), multiply_impl);
  tpool->wait();
}

inline double SparseMatrixChunked::compute_scaling_factor_for_scale(
    const std::vector<double>& weights) const {
  if (empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double sum = 0.0;
  double norm_sum = 0.0;

  for (const auto& matrix : _chunks) {
    for (std::size_t i = 0; i < matrix.size(); ++i) {
      const auto bin1_id = matrix.bin1_ids()[i];
      const auto bin2_id = matrix.bin2_ids()[i];
      const auto count = matrix.counts()[i];

      const auto w1 = weights[bin1_id];
      const auto w2 = weights[bin2_id];

      if (!std::isnan(w1) && !std::isnan(w2)) {
        const auto cfx = bin1_id != bin2_id ? 2.0 : 1.0;
        sum += count * cfx;
        norm_sum += (count * cfx) / (w1 * w2);
      }
    }
  }

  return std::sqrt(norm_sum / sum);
}

inline FileBackedSparseMatrix::FileBackedSparseMatrix(std::filesystem::path tmp_file,
                                                      std::size_t chunk_size, int compression_lvl)
    : _path(std::move(tmp_file)),
      _fs(filestream::FileStream::create(_path.string())),
      _chunk_size(chunk_size),
      _compression_lvl(compression_lvl),
      _zstd_cctx(ZSTD_createCCtx()),
      _zstd_dctx(ZSTD_createDCtx()) {}

inline FileBackedSparseMatrix::~FileBackedSparseMatrix() noexcept {
  try {
    if (!_path.empty() && std::filesystem::exists(_path)) {
      _fs = filestream::FileStream{};
      std::filesystem::remove(_path);
    }
    // NOLINTNEXTLINE
  } catch (...) {
  }
}

inline bool FileBackedSparseMatrix::empty() const noexcept { return size() == 0; }
inline std::size_t FileBackedSparseMatrix::size() const noexcept { return _size; }
inline void FileBackedSparseMatrix::clear(bool shrink_to_fit_) {
  _index.clear();
  _fs = filestream::FileStream{};
  std::filesystem::remove(_path);
  _path = "";
  _size = 0;
  _matrix.clear(shrink_to_fit_);
}

inline void FileBackedSparseMatrix::push_back(std::uint64_t bin1_id, std::uint64_t bin2_id,
                                              double count, std::size_t bin_offset) {
  if (_matrix.size() == _chunk_size) {
    write_chunk();
  }

  _matrix.push_back(bin1_id, bin2_id, count, bin_offset);
  ++_size;
}

inline void FileBackedSparseMatrix::finalize() {
  if (!_matrix.empty()) {
    write_chunk();
  }
  _fs = filestream::FileStream(_path.string());
}

inline void FileBackedSparseMatrix::marginalize(VectorOfAtomicDecimals& marg,
                                                BS::thread_pool* tpool, bool init_buffer) const {
  auto marginalize_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    filestream::FileStream fs(_path.string());
    auto matrix = _matrix;
    std::string buff{};
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(static_cast<std::streamoff>(offset));
      matrix.deserialize(fs, buff, *zstd_dctx);
      matrix.marginalize(marg, false);
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

inline void FileBackedSparseMatrix::marginalize_nnz(VectorOfAtomicDecimals& marg,
                                                    BS::thread_pool* tpool,
                                                    bool init_buffer) const {
  auto marginalize_nnz_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    filestream::FileStream fs(_path.string());
    auto matrix = _matrix;
    std::string buff{};
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(static_cast<std::streamoff>(offset));
      matrix.deserialize(fs, buff, *zstd_dctx);
      matrix.marginalize_nnz(marg, false);
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

inline void FileBackedSparseMatrix::times_outer_product_marg(VectorOfAtomicDecimals& marg,
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
      matrix.times_outer_product_marg(marg, biases, weights, false);
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

inline void FileBackedSparseMatrix::multiply(VectorOfAtomicDecimals& buffer,
                                             nonstd::span<const double> cfx, BS::thread_pool* tpool,
                                             bool init_buffer) const {
  auto multiply_impl = [&](std::size_t istart, std::size_t iend) {
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx(ZSTD_createDCtx());
    filestream::FileStream fs(_path.string());
    auto matrix = _matrix;
    std::string buff{};
    for (const auto offset : nonstd::span(_index).subspan(istart, iend - istart)) {
      fs.seekg(static_cast<std::streamoff>(offset));
      matrix.deserialize(fs, buff, *zstd_dctx);
      matrix.multiply(buffer, cfx, false);
    }
  };

  buffer.resize(cfx.size());
  if (init_buffer) {
    buffer.fill(0);
  }

  if (_index.size() == 1 || !tpool) {
    multiply_impl(0, _index.size());
    return;
  }

  const auto offsets = compute_chunk_offsets(_index.size(), tpool->get_thread_count());

  for (std::size_t i = 1; i < offsets.size(); ++i) {
    const auto i0 = offsets[i - 1];
    const auto i1 = offsets[i];

    tpool->detach_task([&, ii = i0, jj = i1]() { multiply_impl(ii, jj); });
  }
  tpool->wait();
}

inline double FileBackedSparseMatrix::compute_scaling_factor_for_scale(
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

      if (!std::isnan(w1) && !std::isnan(w2)) {
        const auto cfx = bin1_id != bin2_id ? 2.0 : 1.0;
        sum += count * cfx;
        norm_sum += (count * cfx) / (w1 * w2);
      }
    }
  }

  return std::sqrt(norm_sum / sum);
}

inline void FileBackedSparseMatrix::write_chunk() {
  assert(!_matrix.empty());
  _index.push_back(_fs.tellp());
  _matrix.finalize();
  _matrix.serialize(_fs, _buff, *_zstd_cctx, _compression_lvl);
  _matrix.clear();
}

inline std::vector<std::size_t> FileBackedSparseMatrix::compute_chunk_offsets(
    std::size_t size, std::size_t num_chunks) {
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

}  // namespace hictk::balancing::internal
