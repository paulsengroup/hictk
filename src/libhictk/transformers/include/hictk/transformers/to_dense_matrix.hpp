// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_EIGEN

#include <Eigen/Dense>
#include <cstdint>
#include <string_view>
#include <type_traits>

#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/file.hpp"
#include "hictk/hic/pixel_selector.hpp"
#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
class ToDenseMatrix {
 private:
  using PixelIt = decltype(std::declval<PixelSelector>().template begin<N>());
  using PixelT = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  static_assert(std::is_same_v<PixelT, hictk::ThinPixel<N>>);

  PixelSelector _sel{};
  QuerySpan _span{QuerySpan::full};

 public:
  using MatrixT = Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  ToDenseMatrix(PixelSelector&& selector, N n, QuerySpan span = QuerySpan::full);
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

  [[nodiscard]] auto init_matrix() const -> MatrixT;

  [[nodiscard]] std::pair<Eigen::Vector<N, Eigen::Dynamic>, Eigen::Vector<N, Eigen::Dynamic>>
  slice_weights(const cooler::PixelSelector& sel) const;
  [[nodiscard]] std::pair<Eigen::Vector<N, Eigen::Dynamic>, Eigen::Vector<N, Eigen::Dynamic>>
  slice_weights(const hic::PixelSelector& sel) const;
  [[nodiscard]] std::pair<Eigen::Vector<N, Eigen::Dynamic>, Eigen::Vector<N, Eigen::Dynamic>>
  slice_weights(const hic::PixelSelectorAll& sel) const;
  [[nodiscard]] std::pair<Eigen::Vector<N, Eigen::Dynamic>, Eigen::Vector<N, Eigen::Dynamic>>
  slice_weights(const hictk::PixelSelector& sel) const;

  [[nodiscard]] static std::pair<Eigen::Vector<N, Eigen::Dynamic>, Eigen::Vector<N, Eigen::Dynamic>>
  slice_weights(const balancing::Weights& weights1, const balancing::Weights& weights2,
                std::int64_t offset1, std::int64_t offset2, std::int64_t size1, std::int64_t size2);

  void validate_dtype() const;
};

}  // namespace hictk::transformers

#include "./impl/to_dense_matrix_impl.hpp"

#endif
