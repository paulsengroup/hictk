// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <BS_thread_pool.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <nonstd/span.hpp>
#include <utility>
#include <vector>

#include "hictk/balancing/sparse_matrix.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"

namespace hictk::balancing {

class ICE {
  std::vector<std::uint64_t> _chrom_offsets{};
  std::vector<double> _biases{};
  std::vector<double> _variance{};
  std::vector<double> _scale{};

  struct Result {
    double scale{};
    double variance{};
  };

 public:
  enum class Type : std::uint_fast8_t { cis, trans, gw };

  // NOLINTBEGIN(*-avoid-magic-numbers)
  struct Params {
    double tol{1.0e-5};
    std::size_t max_iters{200};
    std::size_t num_masked_diags{2};
    std::size_t min_nnz{10};
    std::size_t min_count{0};
    double mad_max{5.0};
    std::filesystem::path tmpfile{};
    std::size_t chunk_size{10'000'000};
    std::size_t threads{1};
  };

  // NOLINTNEXTLINE(cert-err58-cpp)
  inline static const Params DefaultParams{1.0e-5, 200, 2, 10, 0, 5.0, "", 10'000'000, 1};
  // NOLINTEND(*-avoid-magic-numbers)

  template <typename File>
  explicit ICE(const File& f, Type type = Type::gw, const Params& params = DefaultParams);

  [[nodiscard]] balancing::Weights get_weights(bool rescale = true) const;
  [[nodiscard]] std::vector<double> scale() const noexcept;
  [[nodiscard]] std::vector<double> variance() const noexcept;

 private:
  template <typename File>
  void balance_in_memory(const File& f, Type type, double tol, std::size_t max_iters,
                         std::size_t num_masked_diags, std::size_t min_nnz, std::size_t min_count,
                         double mad_max, BS::light_thread_pool* tpool);

  template <typename File>
  void balance_chunked(const File& f, Type type, double tol, std::size_t max_iters,
                       std::size_t num_masked_diags, std::size_t min_nnz, std::size_t min_count,
                       double mad_max, const std::filesystem::path& tmpfile, std::size_t chunk_size,
                       BS::light_thread_pool* tpool);

  template <typename MatrixT>
  void balance_gw(const MatrixT& matrix, std::size_t max_iters, double tol,
                  BS::light_thread_pool* tpool);

  template <typename MatrixT>
  void balance_cis(const MatrixT& matrix, const Chromosome& chrom, std::size_t max_iters,
                   double tol, BS::light_thread_pool* tpool);

  template <typename MatrixT>
  void balance_trans(const MatrixT& matrix, const BinTable& bins, std::size_t max_iters, double tol,
                     BS::light_thread_pool* tpool);

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix(const File& f, Type type,
                                                    std::size_t num_masked_diags)
      -> internal::SparseMatrixChunked;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_gw(const File& f, std::size_t num_masked_diags)
      -> internal::SparseMatrixChunked;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_cis(const File& f, const Chromosome& chrom,
                                                        std::size_t bin_offset,
                                                        std::size_t num_masked_diags)
      -> internal::SparseMatrixChunked;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_cis(const File& f, std::size_t num_masked_diags)
      -> internal::SparseMatrixChunked;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_trans(const File& f,
                                                          std::size_t num_masked_diags)
      -> internal::SparseMatrixChunked;

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked(const File& f, Type type,
                                                            std::size_t num_masked_diags,
                                                            const std::filesystem::path& tmpfile,
                                                            std::size_t chunk_size)
      -> internal::FileBackedSparseMatrix;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked_gw(const File& f,
                                                               std::size_t num_masked_diags,
                                                               const std::filesystem::path& tmpfile,
                                                               std::size_t chunk_size)
      -> internal::FileBackedSparseMatrix;

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked_cis(
      const File& f, const Chromosome& chrom, std::size_t bin_offset, std::size_t num_masked_diags,
      const std::filesystem::path& tmpfile, std::size_t chunk_size)
      -> internal::FileBackedSparseMatrix;
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked_cis(
      const File& f, std::size_t num_masked_diags, const std::filesystem::path& tmpfile,
      std::size_t chunk_size) -> internal::FileBackedSparseMatrix;

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_chunked_trans(
      const File& f, std::size_t num_masked_diags, const std::filesystem::path& tmpfile,
      std::size_t chunk_size) -> internal::FileBackedSparseMatrix;

  template <typename MatrixT>
  [[nodiscard]] static auto inner_loop(const MatrixT& matrix, nonstd::span<double> biases,
                                       internal::VectorOfAtomicDecimals& marg,
                                       nonstd::span<const double> weights = {},
                                       BS::light_thread_pool* tpool = nullptr) -> Result;
  [[nodiscard]] static std::pair<double, std::size_t> aggregate_marg(
      nonstd::span<const double> marg, BS::light_thread_pool* tpool);

  static void update_biases(nonstd::span<const double> marg, nonstd::span<double> biases,
                            double avg_nzmarg, BS::light_thread_pool* tpool);

  [[nodiscard]] static double compute_ssq_nzmarg(nonstd::span<const double> marg, double avg_nzmarg,
                                                 BS::light_thread_pool* tpool);

  template <typename MatrixT>
  static void min_nnz_filtering(internal::VectorOfAtomicDecimals& marg, const MatrixT& matrix,
                                nonstd::span<double> biases, std::size_t min_nnz,
                                BS::light_thread_pool* tpool);

  static void min_count_filtering(nonstd::span<double> biases, std::size_t min_count,
                                  nonstd::span<const double> marg);

  static void mad_max_filtering(nonstd::span<const std::uint64_t> chrom_offsets,
                                nonstd::span<double> biases, nonstd::span<double> marg,
                                double mad_max);

  template <typename MatrixT>
  static void initialize_biases(const MatrixT& matrix, nonstd::span<double> biases,
                                nonstd::span<const std::uint64_t> chrom_bin_offsets,
                                std::size_t min_nnz, std::size_t min_count, double mad_max,
                                BS::light_thread_pool* tpool);

  [[nodiscard]] static std::vector<double> compute_weights_from_chromosome_sizes(
      const BinTable& bins, nonstd::span<std::uint64_t> chrom_bin_offsets);

  template <typename Vector>
  [[nodiscard]] static bool process_in_parallel(const Vector& marg,
                                                const BS::light_thread_pool* tpool) noexcept;
};

}  // namespace hictk::balancing

#include "./impl/ice_impl.hpp"  //NOLINT
