// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <nonstd/span.hpp>
#include <vector>

#include "hictk/balancing/sparse_matrix.hpp"
#include "hictk/bin_table.hpp"

namespace hictk::balancing {

class ICE {
  std::vector<std::size_t> _chrom_offsets{};
  std::vector<double> _biases{};
  std::vector<double> _variance{};
  std::vector<double> _scale{};

  struct Result {
    double scale{};
    double variance{};
  };

 public:
  enum Type { cis, trans, gw };

  struct Params {
    double tol{1.0e-6};
    std::size_t max_iters{200};
    std::size_t num_masked_diags{2};
    std::size_t min_nnz{10};
    std::size_t min_count{0};
    double mad_max{5.0};
    std::filesystem::path tmpfile{};
    std::size_t chunk_size{10'000'000};
  };

  inline static const Params DefaultParams{1.0e-6, 200, 2, 10, 0, 5.0, "", 10'000'000};  // NOLINT

  template <typename File>
  explicit ICE(const File& f, Type type = Type::gw, const Params& params = DefaultParams);

  [[nodiscard]] std::vector<double> get_weights(bool rescale = true) const;
  [[nodiscard]] std::vector<double> scale() const noexcept;
  [[nodiscard]] std::vector<double> variance() const noexcept;

 private:
  template <typename File>
  void balance_in_memory(const File& f, Type type, double tol, std::size_t max_iters,
                         std::size_t num_masked_diags, std::size_t min_nnz, std::size_t min_count,
                         double mad_max);

  template <typename File>
  void balance_chunked(const File& f, Type type, double tol, std::size_t max_iters,
                       std::size_t num_masked_diags, std::size_t min_nnz, std::size_t min_count,
                       double mad_max, const std::filesystem::path& tmpfile,
                       std::size_t chunk_size);

  template <typename MatrixT>
  void balance_gw(const MatrixT& matrix, std::size_t max_iters, double tol);

  template <typename MatrixT>
  void balance_cis(const MatrixT& matrix, const BinTable& bins, std::size_t max_iters, double tol);

  template <typename MatrixT>
  void balance_trans(const MatrixT& matrix, const BinTable& bins, std::size_t max_iters,
                     double tol);

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix(const File& f, Type type,
                                                    std::size_t num_masked_diags) -> SparseMatrix;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_gw(const File& f, std::size_t num_masked_diags)
      -> SparseMatrix;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_cis(const File& f, std::size_t num_masked_diags)
      -> SparseMatrix;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_trans(const File& f,
                                                          std::size_t num_masked_diags)
      -> SparseMatrix;

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked(const File& f, Type type,
                                                            std::size_t num_masked_diags,
                                                            const std::filesystem::path& tmpfile,
                                                            std::size_t chunk_size)
      -> SparseMatrixChunked;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked_gw(const File& f,
                                                               std::size_t num_masked_diags,
                                                               const std::filesystem::path& tmpfile,
                                                               std::size_t chunk_size)
      -> SparseMatrixChunked;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked_cis(
      const File& f, std::size_t num_masked_diags, const std::filesystem::path& tmpfile,
      std::size_t chunk_size) -> SparseMatrixChunked;

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked_trans(
      const File& f, std::size_t num_masked_diags, const std::filesystem::path& tmpfile,
      std::size_t chunk_size) -> SparseMatrixChunked;

  template <typename MatrixT>
  [[nodiscard]] static auto inner_loop(const MatrixT& matrix, nonstd::span<double> biases,
                                       nonstd::span<const double> weights = {}) -> Result;

  template <typename MatrixT>
  static void min_nnz_filtering(const MatrixT& matrix, nonstd::span<double> biases,
                                std::size_t min_nnz);

  static void min_count_filtering(nonstd::span<double> biases, std::size_t min_count,
                                  nonstd::span<const double> marg);

  static void mad_max_filtering(nonstd::span<const std::size_t> chrom_offsets,
                                nonstd::span<double> biases, nonstd::span<double> marg,
                                double mad_max);

  template <typename MatrixT>
  static void initialize_biases(const MatrixT& matrix, nonstd::span<double> biases,
                                nonstd::span<const std::size_t> chrom_bin_offsets,
                                std::size_t min_nnz, std::size_t min_count, double mad_max);

  [[nodiscard]] static std::vector<std::size_t> read_chrom_bin_offsets(const BinTable& bins);

  [[nodiscard]] static std::vector<double> compute_weights_from_chromosome_sizes(
      const BinTable& bins, nonstd::span<std::size_t> chrom_bin_offsets);
};

}  // namespace hictk::balancing

#include "./impl/ice_mem_impl.hpp"
