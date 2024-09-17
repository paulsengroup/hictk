// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_EIGEN

#include <Eigen/SparseCore>
#include <cstdint>
#include <type_traits>

#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
class ToSparseMatrix {
  using PixelIt = decltype(std::declval<PixelSelector>().template begin<N>());
  using PixelT = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  static_assert(std::is_same_v<PixelT, ThinPixel<N>>);

  PixelSelector _sel{};
  QuerySpan _span{QuerySpan::upper_triangle};

 public:
  using MatrixT = Eigen::SparseMatrix<N, Eigen::RowMajor>;
  ToSparseMatrix(PixelSelector&& selector, N n, QuerySpan span = QuerySpan::upper_triangle);
  [[nodiscard]] auto operator()() -> MatrixT;

 private:
  [[nodiscard]] std::string_view chrom1() const noexcept;
  [[nodiscard]] std::string_view chrom2() const noexcept;

  [[nodiscard]] static std::int64_t num_bins(const PixelCoordinates& coords,
                                             const BinTable& bins) noexcept;
  [[nodiscard]] std::int64_t num_rows() const noexcept;
  [[nodiscard]] std::int64_t num_cols() const noexcept;

  [[nodiscard]] static std::int64_t offset(const PixelCoordinates& coords) noexcept;
  [[nodiscard]] std::int64_t row_offset() const noexcept;
  [[nodiscard]] std::int64_t col_offset() const noexcept;

  void validate_dtype() const;
};

}  // namespace hictk::transformers

#include "./impl/to_sparse_matrix_impl.hpp"

#endif
