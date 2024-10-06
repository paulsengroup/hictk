// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_EIGEN

#include <Eigen/SparseCore>
#include <cstdint>
#include <memory>
#include <type_traits>
#include <vector>

#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
class ToSparseMatrix {
  using PixelIt = decltype(std::declval<PixelSelector>().template begin<N>());
  using PixelT = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  static_assert(std::is_same_v<PixelT, ThinPixel<N>>);

  using MatrixRowMajor = Eigen::SparseMatrix<N, Eigen::RowMajor>;
  using MatrixColMajor = Eigen::SparseMatrix<N, Eigen::ColMajor>;

  std::shared_ptr<const PixelSelector> _sel{};
  QuerySpan _span{QuerySpan::upper_triangle};

 public:
  using MatrixT = MatrixRowMajor;
  ToSparseMatrix() = delete;
  ToSparseMatrix(PixelSelector selector, N n, QuerySpan span = QuerySpan::upper_triangle);
  ToSparseMatrix(std::shared_ptr<const PixelSelector> selector, N n,
                 QuerySpan span = QuerySpan::upper_triangle);

  ToSparseMatrix(const ToSparseMatrix& other) = delete;
  ToSparseMatrix(ToSparseMatrix&& other) noexcept = default;

  ~ToSparseMatrix() noexcept = default;

  ToSparseMatrix& operator=(const ToSparseMatrix& other) = delete;
  ToSparseMatrix& operator=(ToSparseMatrix&& other) noexcept = default;

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

  template <typename Matrix, typename SetterOp>
  [[nodiscard]] static auto fill_row(PixelIt first_pixel, PixelIt last_pixel,
                                     MatrixRowMajor& matrix_ut, Matrix& matrix_lt,
                                     std::vector<std::int64_t>& row_sizes,
                                     std::vector<ThinPixel<N>>& buffer, bool symmetric_upper,
                                     std::int64_t offset1, std::int64_t offset2,
                                     bool populate_lower_triangle, bool populate_upper_triangle,
                                     SetterOp matrix_setter) -> PixelIt;
};

}  // namespace hictk::transformers

#include "./impl/to_sparse_matrix_impl.hpp"

#endif
