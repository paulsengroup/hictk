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
  std::vector<std::uint64_t> _chrom_offsets{};
  std::vector<double> _biases{};
  std::vector<double> _variance{};
  std::vector<double> _scale{};

  struct Result {
    std::vector<std::uint64_t> offsets{};
    std::vector<double> scales{};
    std::vector<double> weights{};
  };

 public:
  enum Type { cis, trans, gw };

  struct Params {
    double tol{1.0e-4};
    std::size_t max_iters{200};
    double max_percentile{10};
    double erez{1.0e-5};
    double max_row_sum_error = 0.05;
    double delta = 0.05;
  };

  // NOLINTNEXTLINE
  inline static const Params DefaultParams{1.0e-4, 200, 10.0, 1.0e-5, 0.05, 0.05};

  template <typename File>
  explicit SCALE(const File& f, Type type = Type::gw, const Params& params = DefaultParams);
  template <typename PixelIt>
  SCALE(PixelIt first, PixelIt last, const BinTable& bins, const Params& params = DefaultParams);

  [[nodiscard]] std::vector<double> get_weights(bool rescale = true) const;
  [[nodiscard]] const std::vector<double>& get_scale() const noexcept;

 private:
  template <typename File>
  [[nodiscard]] static auto compute_cis(const File& f, const Params& params) -> Result;
  template <typename File>
  [[nodiscard]] static auto compute_trans(const File& f, const Params& params) -> Result;
  template <typename File>
  [[nodiscard]] static auto compute_gw(const File& f, const Params& params) -> Result;

  [[nodiscard]] static VC::Type map_type_to_vc(Type type) noexcept;

  [[nodiscard]] static std::vector<double> matrix_vect_mult(
      const std::vector<ThinPixel<double>>& pixels, const std::vector<double>& cfx,
      std::size_t offset, std::size_t i0 = 0, std::size_t i1 = 0);
  static void matrix_vect_mult(const std::vector<ThinPixel<double>>& pixels,
                               const std::vector<double>& cfx, std::vector<double>& sum_vect,
                               std::size_t offset, std::size_t i0 = 0, std::size_t i1 = 0) noexcept;
  static std::vector<double> update_weights(const std::vector<bool>& bad,
                                            std::vector<double>& weights,
                                            const std::vector<double>& target,
                                            std::vector<double>& d_vector,
                                            const std::vector<ThinPixel<double>>& pixels,
                                            std::size_t offset) noexcept;

  static void geometric_mean(const std::vector<double>& v1, const std::vector<double>& v2,
                             std::vector<double>& vout) noexcept;

  [[nodiscard]] static std::pair<double, std::uint64_t> compute_convergence_error(
      const std::vector<double>& calculated_vector_b, const std::vector<double>& current,
      const std::vector<bool>& bad, double tolerance) noexcept;

  [[nodiscard]] static double compute_final_error(const std::vector<double>& col,
                                                  const std::vector<double>& scale,
                                                  const std::vector<double>& target,
                                                  const std::vector<bool>& bad) noexcept;
  static void multiply(std::vector<double>& v1, const std::vector<double>& v2) noexcept;
  [[nodiscard]] static double compute_scale(const std::vector<ThinPixel<double>>& pixels,
                                            const std::vector<double>& weights,
                                            std::size_t offset) noexcept;
};

}  // namespace hictk::balancing

#include "./impl/scale_impl.hpp"  //NOLINT
