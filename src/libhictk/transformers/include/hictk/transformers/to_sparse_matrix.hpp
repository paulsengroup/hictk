// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_EIGEN

#include <Eigen/SparseCore>
#include <cstdint>
#include <type_traits>

#include "hictk/type_traits.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
class ToSparseMatrix {
  using PixelIt = decltype(std::declval<PixelSelector>().template begin<N>());
  using PixelT = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  static_assert(std::is_same_v<PixelT, hictk::ThinPixel<N>>);

  PixelSelector _sel{};

 public:
  ToSparseMatrix(PixelSelector&& selector, N n);
  [[nodiscard]] auto operator()() -> Eigen::SparseMatrix<N>;

 private:
  [[nodiscard]] std::int64_t num_rows() const noexcept;
  [[nodiscard]] std::int64_t num_cols() const noexcept;

  [[nodiscard]] std::uint64_t row_offset() const noexcept;
  [[nodiscard]] std::uint64_t col_offset() const noexcept;
};

}  // namespace hictk::transformers

#include "./impl/to_sparse_matrix_impl.hpp"

#endif
