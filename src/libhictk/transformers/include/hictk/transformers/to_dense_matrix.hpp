// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_EIGEN

#include <Eigen/Dense>
#include <cstdint>
#include <string_view>
#include <type_traits>

#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
class ToDenseMatrix {
  using PixelIt = decltype(std::declval<PixelSelector>().template begin<N>());
  using PixelT = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  static_assert(std::is_same_v<PixelT, hictk::ThinPixel<N>>);

  PixelSelector _sel{};
  bool _mirror{};

 public:
  ToDenseMatrix(PixelSelector&& selector, N n, bool mirror = true);
  [[nodiscard]] auto operator()() -> Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic>;

 private:
  [[nodiscard]] std::string_view chrom1() const noexcept;
  [[nodiscard]] std::string_view chrom2() const noexcept;
  [[nodiscard]] std::int64_t num_rows() const noexcept;
  [[nodiscard]] std::int64_t num_cols() const noexcept;

  [[nodiscard]] std::uint64_t row_offset() const noexcept;
  [[nodiscard]] std::uint64_t col_offset() const noexcept;
};

}  // namespace hictk::transformers

#include "./impl/to_dense_matrix_impl.hpp"

#endif
