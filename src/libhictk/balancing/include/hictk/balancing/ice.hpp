// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <variant>
#include <vector>

namespace hictk::balancing {

class ICE {
  std::vector<double> _biases{};
  double _variance{0.0};
  double _scale{std::numeric_limits<double>::quiet_NaN()};
  std::variant<std::int64_t, double> _sum{};

  struct Result {
    double scale;
    double variance;
  };

 public:
  template <typename PixelIt>
  ICE(PixelIt first_pixel, PixelIt last_pixel, std::size_t num_rows, double tol = 1.0e-5,
      std::size_t max_iters = 200, std::size_t num_masked_diags = 2, std::size_t min_nnz = 10,
      double min_count = 0);

  [[nodiscard]] std::vector<double> get_weights(bool rescale = true) const;
  [[nodiscard]] double scale() const noexcept;
  [[nodiscard]] double variance() const noexcept;

 private:
  template <typename PixelIt>
  static std::tuple<std::vector<std::size_t>, std::vector<std::size_t>, std::vector<double>>
  construct_sparse_matrix(PixelIt first_pixel, PixelIt last_pixel, std::size_t num_masked_diags);

  [[nodiscard]] static auto inner_loop(const std::vector<std::size_t>& bin1_ids,
                                       const std::vector<std::size_t>& bin2_ids,
                                       std::vector<double> counts, std::vector<double>& biases,
                                       std::vector<double>& marg_buffer) -> Result;

  static void times_outer_product(const std::vector<std::size_t>& bin1_ids,
                                  const std::vector<std::size_t>& bin2_ids,
                                  std::vector<double>& counts, const std::vector<double>& biases);

  static void marginalize(const std::vector<std::size_t>& bin1_ids,
                          const std::vector<std::size_t>& bin2_ids, std::vector<double>& counts,
                          std::vector<double>& marg);

  static void filter_rows_by_nnz(const std::vector<std::size_t>& bin1_ids,
                                 const std::vector<std::size_t>& bin2_ids,
                                 std::vector<double> counts, std::vector<double>& biases,
                                 std::size_t min_nnz, std::vector<double>& marg_buff);

};

}  // namespace hictk::balancing

#include "../../../ice_impl.hpp"
