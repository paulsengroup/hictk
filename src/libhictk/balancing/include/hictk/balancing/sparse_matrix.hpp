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

#include "hictk/binary_buffer.hpp"
#include "hictk/common.hpp"
#include "hictk/default_delete.hpp"
#include "hictk/filestream.hpp"

namespace hictk::balancing::internal {

class AtomicBitSet {
  using I = std::uint8_t;
  std::vector<std::atomic<I>> _buff{};
  std::size_t _size{};

 public:
  AtomicBitSet() = default;
  explicit AtomicBitSet(std::size_t size_, bool value = false);
  AtomicBitSet(const AtomicBitSet& other);
  AtomicBitSet(AtomicBitSet&& other) = default;

  ~AtomicBitSet() = default;

  AtomicBitSet& operator=(const AtomicBitSet& other);
  AtomicBitSet& operator=(AtomicBitSet&& other) = default;

  void atomic_set(std::size_t i, bool value) noexcept;
  [[nodiscard]] bool atomic_test(std::size_t i) const noexcept;

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;

  void fill(bool value) noexcept;
  void resize(std::size_t size_, bool value = false);

 private:
  [[nodiscard]] static std::size_t compute_offset(std::size_t i) noexcept;
  [[nodiscard]] static bool atomic_test(const std::vector<std::atomic<I>>& buff,
                                        std::size_t i) noexcept;
};

class VectorOfAtomicDecimals {
  using I = std::uint64_t;
  using N = std::atomic<I>;
  std::vector<N> _margsi{};
  mutable std::vector<double> _margsd{};
  AtomicBitSet _nanmask{};
  AtomicBitSet _infmask{};
  std::uint64_t _cfxi{};
  double _cfxd{};

  static constexpr std::uint8_t DEFAULT_DECIMAL_BITS = 30;
  double _max_value{compute_max_value(DEFAULT_DECIMAL_BITS)};

 public:
  VectorOfAtomicDecimals() = default;
  explicit VectorOfAtomicDecimals(std::size_t size_ = 0,
                                  std::uint64_t decimal_bits = DEFAULT_DECIMAL_BITS);

  VectorOfAtomicDecimals(const VectorOfAtomicDecimals& other);
  VectorOfAtomicDecimals(VectorOfAtomicDecimals&& other) noexcept = default;

  ~VectorOfAtomicDecimals() = default;

  VectorOfAtomicDecimals& operator=(const VectorOfAtomicDecimals& other);
  VectorOfAtomicDecimals& operator=(VectorOfAtomicDecimals&& other) noexcept = default;

  [[nodiscard]] double operator[](std::size_t i) const noexcept;
  void atomic_add(std::size_t i, double n) noexcept;
  void set(std::size_t i, double n) noexcept;
  void multiply(const std::vector<double>& v) noexcept;

  [[nodiscard]] const std::vector<double>& operator()() const noexcept;
  [[nodiscard]] std::vector<double>& operator()() noexcept;

  void fill(double value = 0) noexcept;
  void resize(std::size_t size_, double value = 0);

  [[nodiscard]] std::uint8_t decimal_bits() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::pair<double, double> domain(bool include_inf = true) const noexcept;

 private:
  [[nodiscard]] auto encode(double n) const noexcept -> I;
  [[nodiscard]] double decode(I n) const noexcept;
  [[nodiscard]] constexpr bool overflows(double n) const noexcept;
  [[nodiscard]] static double compute_max_value(std::uint8_t decimal_bits) noexcept;
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
  void reserve(std::size_t capacity);
  void shrink_to_fit() noexcept;
  void finalize();

  [[nodiscard]] const std::vector<std::uint64_t>& bin1_ids() const noexcept;
  [[nodiscard]] const std::vector<std::uint64_t>& bin2_ids() const noexcept;
  [[nodiscard]] const std::vector<double>& counts() const noexcept;

  void push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count,
                 std::size_t bin_offset = 0);

  void serialize(filestream::FileStream<>& fs, std::string& tmpbuff, ZSTD_CCtx& ctx,
                 int compression_lvl = 3) const;
  void deserialize(filestream::FileStream<>& fs, std::string& tmpbuff, ZSTD_DCtx& ctx);

  void marginalize(VectorOfAtomicDecimals& marg, bool init_buffer = true) const;
  void marginalize_nnz(VectorOfAtomicDecimals& marg, bool init_buffer = true) const;
  void times_outer_product_marg(VectorOfAtomicDecimals& marg, nonstd::span<const double> biases,
                                nonstd::span<const double> weights, bool init_buffer = true) const;

  void multiply(VectorOfAtomicDecimals& buffer, nonstd::span<const double> cfx,
                bool init_buffer = true) const;

  [[nodiscard]] double compute_scaling_factor_for_scale(const std::vector<double>& weights) const;
};

class SparseMatrixChunked {
  std::vector<SparseMatrix> _chunks{};
  std::size_t _size{};
  std::size_t _chunk_size{};

 public:
  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  explicit SparseMatrixChunked(std::size_t chunk_size = 16UL << 20U);

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::size_t num_chunks() const noexcept;
  [[nodiscard]] std::size_t chunk_size() const noexcept;
  void shrink_to_fit() noexcept;
  void clear(bool shrink_to_fit_ = false);

  void push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count,
                 std::size_t bin_offset = 0);
  void finalize();

  void marginalize(VectorOfAtomicDecimals& marg, BS::light_thread_pool* tpool = nullptr,
                   bool init_buffer = true) const;
  void marginalize_nnz(VectorOfAtomicDecimals& marg, BS::light_thread_pool* tpool = nullptr,
                       bool init_buffer = true) const;
  void times_outer_product_marg(VectorOfAtomicDecimals& marg, nonstd::span<const double> biases,
                                nonstd::span<const double> weights,
                                BS::light_thread_pool* tpool = nullptr,
                                bool init_buffer = true) const;

  void multiply(VectorOfAtomicDecimals& buffer, nonstd::span<const double> cfx,
                BS::light_thread_pool* tpool = nullptr, bool init_buffer = true) const;

  [[nodiscard]] double compute_scaling_factor_for_scale(const std::vector<double>& weights) const;
};

class FileBackedSparseMatrix {
  mutable SparseMatrix _matrix{};
  mutable std::string _buff{};
  std::filesystem::path _path{};
  mutable filestream::FileStream<> _fs{};

  std::vector<std::streampos> _index{};
  std::size_t _size{};
  std::size_t _chunk_size{};
  int _compression_lvl{};

  std::unique_ptr<ZSTD_CCtx_s> _zstd_cctx{};
  std::unique_ptr<ZSTD_DCtx_s> _zstd_dctx{};

 public:
  FileBackedSparseMatrix() = default;
  FileBackedSparseMatrix(std::filesystem::path tmp_file, std::size_t chunk_size,
                         int compression_lvl = 3);

  FileBackedSparseMatrix(const FileBackedSparseMatrix& other) = delete;
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ < 10
  FileBackedSparseMatrix(FileBackedSparseMatrix&& other) = default;
#elif defined(__clang__) && __clang__ < 9
  FileBackedSparseMatrix(FileBackedSparseMatrix&& other) = default;
#else
  FileBackedSparseMatrix(FileBackedSparseMatrix&& other) noexcept = default;
#endif

  ~FileBackedSparseMatrix() noexcept;

  FileBackedSparseMatrix& operator=(const FileBackedSparseMatrix& other) = delete;
#if defined(__GNUC__) && defined(__clang__) && __clang_major__ > 8
  FileBackedSparseMatrix& operator=(FileBackedSparseMatrix&& other) noexcept = default;
#elif defined(__GNUC__) && __GNUC__ > 9
  FileBackedSparseMatrix& operator=(FileBackedSparseMatrix&& other) noexcept = default;
#else
  FileBackedSparseMatrix& operator=(FileBackedSparseMatrix&& other) = default;
#endif

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  void clear(bool shrink_to_fit_ = false);

  void push_back(std::uint64_t bin1_id, std::uint64_t bin2_id, double count,
                 std::size_t bin_offset = 0);
  void finalize();

  void marginalize(VectorOfAtomicDecimals& marg, BS::light_thread_pool* tpool = nullptr,
                   bool init_buffer = true) const;
  void marginalize_nnz(VectorOfAtomicDecimals& marg, BS::light_thread_pool* tpool = nullptr,
                       bool init_buffer = true) const;
  void times_outer_product_marg(VectorOfAtomicDecimals& marg, nonstd::span<const double> biases,
                                nonstd::span<const double> weights,
                                BS::light_thread_pool* tpool = nullptr,
                                bool init_buffer = true) const;

  void multiply(VectorOfAtomicDecimals& buffer, nonstd::span<const double> cfx,
                BS::light_thread_pool* tpool = nullptr, bool init_buffer = true) const;

  [[nodiscard]] double compute_scaling_factor_for_scale(const std::vector<double>& weights) const;

 private:
  void write_chunk();
  [[nodiscard]] static std::vector<std::size_t> compute_chunk_offsets(std::size_t size,
                                                                      std::size_t num_chunks);
};

}  // namespace hictk::balancing::internal

#include "./impl/sparse_matrix_impl.hpp"  // NOLINT
