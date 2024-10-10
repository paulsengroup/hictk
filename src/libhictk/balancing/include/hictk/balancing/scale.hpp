// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
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
#include "hictk/balancing/vc.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"

namespace hictk::balancing {

class SCALE {
  struct Result {
    std::vector<std::uint64_t> offsets{};
    std::vector<double> scales{};
    std::vector<double> weights{};
  };

  // NOLINTBEGIN(*-avoid-magic-numbers)
  struct ConvergenceStats {
    bool converged{false};
    bool diverged{false};
    std::uint32_t low_convergence{1000};
    std::uint32_t low_divergence{0};
    double error{10.0};
  };
  // NOLINTEND(*-avoid-magic-numbers)

  enum class ControlFlow : std::uint_fast8_t { break_loop, continue_loop };

  // buffers to store final results
  std::vector<std::uint64_t> _chrom_offsets{};
  std::vector<double> _biases{};
  std::vector<double> _variance{};
  std::vector<double> _scale{};

  // intermediate buffers
  bool _yes{false};
  ConvergenceStats _convergence_stats{};
  std::queue<double> _error_queue_iter{};

  std::vector<bool> _bad{};
  std::vector<double> _one{};
  std::vector<double> _biases1{};  // calculated_vector_b
  std::vector<double> _z_target_vector{};
  std::vector<std::uint64_t> _row_wise_nnz{};
  std::uint64_t _nnz_rows{};
  std::uint64_t _low_cutoff{};
  std::uint64_t _upper_bound{};

  std::vector<bool> _bad_conv{};
  std::vector<double> _b_conv{};
  double _ber_conv{};

  std::size_t _iter{};
  std::size_t _tot_iter{};
  std::size_t _max_tot_iters{};

  std::unique_ptr<BS::thread_pool> _tpool{};

 public:
  enum class Type : std::uint_fast8_t { cis, trans, gw };

  // NOLINTBEGIN(*-avoid-magic-numbers)
  struct Params {
    double tol{1.0e-4};
    std::size_t max_iters{200};
    double max_percentile{10};
    double frac_bad_cutoff{1.0e-5};
    double max_row_sum_error{0.05};
    double delta{0.05};
    std::filesystem::path tmpfile{};
    std::size_t chunk_size{10'000'000};
    std::size_t threads{1};
  };

  // NOLINTNEXTLINE(cert-err58-cpp)
  inline static const Params DefaultParams{1.0e-4, 200, 10.0,       1.0e-5, 0.05,
                                           0.05,   "",  10'000'000, 1};
  // NOLINTEND(*-avoid-magic-numbers)

  template <typename File>
  explicit SCALE(const File& f, Type type = Type::gw, const Params& params = DefaultParams);
  template <typename PixelIt>
  SCALE(PixelIt first, PixelIt last, const BinTable& bins, const Params& params = DefaultParams);

  [[nodiscard]] Weights get_weights(bool rescale = true) const;
  [[nodiscard]] const std::vector<double>& get_scale() const noexcept;

 private:
  template <typename File>
  [[nodiscard]] static auto compute_cis(const File& f, const Params& params) -> Result;
  template <typename File>
  [[nodiscard]] static auto compute_trans(const File& f, const Params& params) -> Result;
  template <typename File>
  [[nodiscard]] static auto compute_gw(const File& f, const Params& params) -> Result;

  template <typename Matrix>
  void balance(const Matrix& m, const BinTable& bins, const Params& params);

  [[nodiscard]] static VC::Type map_type_to_vc(Type type) noexcept;

  template <typename Matrix>
  static void update_weights(internal::VectorOfAtomicDecimals& buffer, const std::vector<bool>& bad,
                             internal::VectorOfAtomicDecimals& weights,
                             const std::vector<double>& target, std::vector<double>& d_vector,
                             const Matrix& m, BS::thread_pool* tpool);

  static void geometric_mean(const std::vector<double>& v1, const std::vector<double>& v2,
                             std::vector<double>& vout) noexcept;

  [[nodiscard]] static std::pair<double, std::uint64_t> compute_convergence_error(
      const std::vector<double>& biases, const std::vector<double>& current,
      const std::vector<bool>& bad, double tolerance) noexcept;

  [[nodiscard]] static double compute_final_error(const internal::VectorOfAtomicDecimals& col,
                                                  const std::vector<double>& scale,
                                                  const std::vector<double>& target,
                                                  const std::vector<bool>& bad) noexcept;
  static void multiply(std::vector<double>& v1, const std::vector<double>& v2) noexcept;

  template <typename PixelIt>
  [[nodiscard]] std::variant<internal::SparseMatrixChunked, internal::FileBackedSparseMatrix>
  mask_bins_and_init_buffers(PixelIt first, PixelIt last, std::size_t offset, double max_percentile,
                             const std::filesystem::path& tmpfile, std::size_t chunk_size);
  template <typename Matrix>
  [[nodiscard]] auto handle_convergenece(const Matrix& m, std::vector<double>& dr,
                                         std::vector<double>& dc,
                                         internal::VectorOfAtomicDecimals& row) -> ControlFlow;

  template <typename Matrix>
  [[nodiscard]] auto handle_almost_converged(const Matrix& m, const std::vector<double>& b0,
                                             std::vector<double>& dr, std::vector<double>& dc,
                                             internal::VectorOfAtomicDecimals& row,
                                             double tolerance) -> ControlFlow;

  template <typename Matrix>
  [[nodiscard]] auto handle_diverged(const Matrix& m, const std::vector<double>& b0,
                                     std::vector<double>& dr, std::vector<double>& dc,
                                     internal::VectorOfAtomicDecimals& row, double frac_bad,
                                     double frac_bad_cutoff, double tolerance) -> ControlFlow;

  template <typename PixelIt>
  static std::variant<internal::SparseMatrixChunked, internal::FileBackedSparseMatrix> init_matrix(
      PixelIt first, PixelIt last, std::size_t offset, const std::filesystem::path& tmpfile,
      std::size_t chunk_size);

  [[nodiscard]] std::size_t size() const noexcept;

  void reset_iter() noexcept;
};

}  // namespace hictk::balancing

#include "./impl/scale_impl.hpp"  //NOLINT
