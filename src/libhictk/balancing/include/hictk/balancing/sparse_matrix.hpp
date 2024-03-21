// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <zstd.h>

#include <BS_thread_pool.hpp>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ios>
#include <memory>
#include <nonstd/span.hpp>
#include <string>
#include <type_traits>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/default_delete.hpp"

namespace hictk::balancing {

class VectorOfAtomicDecimals {
  using I = std::uint64_t;
  using N = std::atomic<I>;
  std::vector<N> _margsi{};
  mutable std::vector<double> _margsd{};
  std::uint64_t _cfx{};
  const static auto DEFAULT_DECIMAL_DIGITS = 9ULL;

 public:
  VectorOfAtomicDecimals() = delete;
  explicit VectorOfAtomicDecimals(std::size_t size_ = 0,
                                  std::size_t decimals = DEFAULT_DECIMAL_DIGITS);

  VectorOfAtomicDecimals(const VectorOfAtomicDecimals& other);
  VectorOfAtomicDecimals(VectorOfAtomicDecimals&& other) noexcept = default;

  ~VectorOfAtomicDecimals() = default;

  VectorOfAtomicDecimals& operator=(const VectorOfAtomicDecimals& other);
  VectorOfAtomicDecimals& operator=(VectorOfAtomicDecimals&& other) noexcept = default;

  [[nodiscard]] double operator[](std::size_t i) const noexcept;
  void add(std::size_t i, double n) noexcept;
  void set(std::size_t i, double n) noexcept;
  // multiply is not atomic
  void multiply(const std::vector<double>& v) noexcept;

  [[nodiscard]] const std::vector<double>& operator()() const noexcept;
  [[nodiscard]] std::vector<double>& operator()() noexcept;

  void fill(double value = 0) noexcept;
  void resize(std::size_t size_);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;

 private:
  auto encode(double n) const noexcept -> I;
  double decode(I n) const noexcept;
};

class SparseMatrix {
  std::vector<std::uint64_t> _bin1_ids{};
  std::vector<std::uint64_t> _bin2_ids{};
  std::vector<double> _counts{};

 public:
  SparseMatrix() = default;

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  void clear(bool shrink_to_fit_ = false) noexcept;
  void shrink_to_fit() noexcept;
  void finalize();

  [[nodiscard]] const std::vector<std::uint64_t>& bin1_ids() const noexcept;
  [[nodiscard]] const std::vector<std::uint64_t>& bin2_ids() const noexcept;
  [[nodiscard]] const std::vector<double>& counts() const noexcept;

  void push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count,
                 std::size_t bin_offset = 0);

  void serialize(std::fstream& fs, ZSTD_CCtx& ctx, int compression_lvl = 3) const;
  void deserialize(std::fstream& fs, ZSTD_DCtx& ctx);

  void marginalize(VectorOfAtomicDecimals& marg, BS::thread_pool* tpool = nullptr,
                   bool init_buffer = true) const;
  void marginalize_nnz(VectorOfAtomicDecimals& marg, BS::thread_pool* tpool = nullptr,
                       bool init_buffer = true) const;
  void times_outer_product_marg(VectorOfAtomicDecimals& marg, nonstd::span<const double> biases,
                                nonstd::span<const double> weights,
                                BS::thread_pool* tpool = nullptr, bool init_buffer = true) const;

  void multiply(VectorOfAtomicDecimals& buffer, nonstd::span<const double> cfx,
                BS::thread_pool* tpool = nullptr, bool init_buffer = true) const;

  [[nodiscard]] double compute_scaling_factor_for_scale(const std::vector<double>& weights) const;
};

class SparseMatrixChunked {
  mutable SparseMatrix _matrix{};
  mutable std::string _buff{};
  std::filesystem::path _path{};
  mutable std::fstream _fs{};

  std::vector<std::streamoff> _index{};
  std::size_t _size{};
  std::size_t _chunk_size{};
  int _compression_lvl{};

  std::unique_ptr<ZSTD_CCtx_s> _zstd_cctx{};
  std::unique_ptr<ZSTD_DCtx_s> _zstd_dctx{};

 public:
  SparseMatrixChunked() = default;
  SparseMatrixChunked(std::filesystem::path tmp_file, std::size_t chunk_size,
                      int compression_lvl = 3);

  SparseMatrixChunked(const SparseMatrixChunked& other) = delete;
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 10
  SparseMatrixChunked(SparseMatrixChunked&& other) = default;
#elif defined(__clang__) && __clang__ < 9
  SparseMatrixChunked(SparseMatrixChunked&& other) = default;
#else
  SparseMatrixChunked(SparseMatrixChunked&& other) noexcept = default;
#endif

  ~SparseMatrixChunked() noexcept;

  SparseMatrixChunked& operator=(const SparseMatrixChunked& other) = delete;
#if defined(__GNUC__) && defined(__clang__) && __clang_major__ > 8
  SparseMatrixChunked& operator=(SparseMatrixChunked&& other) noexcept = default;
#elif defined(__GNUC__) && __GNUC__ > 9
  SparseMatrixChunked& operator=(SparseMatrixChunked&& other) noexcept = default;
#else
  SparseMatrixChunked& operator=(SparseMatrixChunked&& other) = default;
#endif

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  void clear(bool shrink_to_fit_ = false);

  void push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count,
                 std::size_t bin_offset = 0);
  void finalize();

  void marginalize(VectorOfAtomicDecimals& marg, BS::thread_pool* tpool = nullptr,
                   bool init_buffer = true) const;
  void marginalize_nnz(VectorOfAtomicDecimals& marg, BS::thread_pool* tpool = nullptr,
                       bool init_buffer = true) const;
  void times_outer_product_marg(VectorOfAtomicDecimals& marg, nonstd::span<const double> biases,
                                nonstd::span<const double> weights,
                                BS::thread_pool* tpool = nullptr, bool init_buffer = true) const;

  void multiply(VectorOfAtomicDecimals& buffer, nonstd::span<const double> cfx,
                BS::thread_pool* tpool = nullptr, bool init_buffer = true) const;

  [[nodiscard]] double compute_scaling_factor_for_scale(const std::vector<double>& weights) const;

 private:
  void write_chunk();
  [[nodiscard]] static std::vector<std::size_t> compute_chunk_offsets(std::size_t size,
                                                                      std::size_t num_chunks);
};

}  // namespace hictk::balancing

#include "./impl/sparse_matrix_impl.hpp"  // NOLINT
