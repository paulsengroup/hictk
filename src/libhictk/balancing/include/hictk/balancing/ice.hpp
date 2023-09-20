// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <nonstd/span.hpp>
#include <vector>

#include "hictk/bin_table.hpp"

namespace hictk::balancing {

struct SparseMatrix {
  std::vector<std::size_t> bin1_ids{};  // NOLINT
  std::vector<std::size_t> bin2_ids{};  // NOLINT
  std::vector<double> counts{};         // NOLINT

  std::vector<std::size_t> chrom_offsets{};  // NOLINT

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  void clear() noexcept;
};

class ICE {
  std::vector<std::size_t> _chrom_bin_offsets{};
  std::vector<double> _biases{};
  std::vector<double> _variance{};
  std::vector<double> _scale{};

  struct Result {
    double scale;
    double variance;
  };

 public:
  enum Type { cis, trans, gw };

  template <typename File>
  ICE(const File& f, Type type = Type::gw, double tol = 1.0e-5, std::size_t max_iters = 200,
      std::size_t num_masked_diags = 2, std::size_t min_nnz = 10, std::size_t min_count = 0,
      double mad_max = 5.0);

  template <typename PixelIt>
  ICE(PixelIt first_pixel, PixelIt last_pixel, const BinTable& bins, double tol = 1.0e-5,
      std::size_t max_iters = 200, std::size_t num_masked_diags = 2, std::size_t min_nnz = 10,
      std::size_t min_count = 0, double mad_max = 5.0);

  [[nodiscard]] std::vector<double> get_weights(bool rescale = true) const;
  [[nodiscard]] std::vector<double> scale() const noexcept;
  [[nodiscard]] std::vector<double> variance() const noexcept;

 private:
  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix(const File& f, Type type,
                                                    std::size_t num_masked_diags,
                                                    std::size_t bin_offset = 0) -> SparseMatrix;

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_gw(const File& f, std::size_t num_masked_diags,
                                                       std::size_t bin_offset) -> SparseMatrix;

  template <typename File>
  [[nodiscard]] static auto construct_sparse_matrix_cis(const File& f, std::size_t num_masked_diags,
                                                        std::size_t bin_offset) -> SparseMatrix;

  [[nodiscard]] static auto inner_loop(nonstd::span<const std::size_t> bin1_ids,
                                       nonstd::span<const std::size_t> bin2_ids,
                                       nonstd::span<const double> counts,
                                       nonstd::span<double> biases,
                                       nonstd::span<double> marg_buffer, std::size_t bin_offset = 0,
                                       nonstd::span<double> weights = {}) -> Result;

  static void marginalize(nonstd::span<const std::size_t> bin1_ids,
                          nonstd::span<const std::size_t> bin2_ids,
                          nonstd::span<const double> counts, nonstd::span<double> marg,
                          std::size_t bin_offset = 0);

  static void times_outer_product_marg(nonstd::span<const std::size_t> bin1_ids,
                                       nonstd::span<const std::size_t> bin2_ids,
                                       nonstd::span<const double> counts,
                                       nonstd::span<const double> biases, nonstd::span<double> marg,
                                       std::size_t bin_offset,
                                       nonstd::span<const double> weights = {});

  static void marginalize_nnz(nonstd::span<const std::size_t> bin1_ids,
                              nonstd::span<const std::size_t> bin2_ids,
                              nonstd::span<const double> counts, nonstd::span<double> marg,
                              std::size_t bin_offset = 0);

  static void min_nnz_filtering(nonstd::span<const std::size_t> bin1_ids,
                                nonstd::span<const std::size_t> bin2_ids,
                                nonstd::span<const double> counts, nonstd::span<double> biases,
                                std::size_t min_nnz, nonstd::span<double> marg_buff,
                                std::size_t bin_offset = 0);
  static void min_count_filtering(nonstd::span<double> biases, std::size_t min_count,
                                  nonstd::span<double> marg);

  static void mad_max_filtering(nonstd::span<const std::size_t> chrom_offsets,
                                nonstd::span<double> biases, std::vector<double> marg,
                                double mad_max);

  static void initialize_biases(nonstd::span<const std::size_t> bin1_ids,
                                nonstd::span<const std::size_t> bin2_ids,
                                nonstd::span<const double> counts, nonstd::span<double> biases,
                                nonstd::span<const std::size_t> chrom_bin_offsets,
                                std::size_t min_nnz, std::size_t min_count, double mad_max);
  [[nodiscard]] static std::vector<std::size_t> read_chrom_bin_offsets(const BinTable& bins);

  [[nodiscard]] static std::vector<double> compute_weights_from_chromosome_sizes(
      const BinTable& bins, nonstd::span<std::size_t> chrom_bin_offsets);

  static void mask_cis_interactions(nonstd::span<const std::size_t> bin1_ids,
                                    nonstd::span<const std::size_t> bin2_ids,
                                    nonstd::span<double> counts,
                                    nonstd::span<const std::size_t> chrom_bin_offsets);
};

}  // namespace hictk::balancing

#include "../../../ice_impl.hpp"
