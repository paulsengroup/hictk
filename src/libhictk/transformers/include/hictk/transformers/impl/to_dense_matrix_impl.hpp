// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#if __has_include(<Eigen/Dense>)
#include <Eigen/Dense>
#else
#include <eigen3/Eigen/Dense>
#endif

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <optional>
#include <utility>
#include <variant>

#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
inline ToDenseMatrix<N, PixelSelector>::ToDenseMatrix(
    PixelSelector sel, [[maybe_unused]] N n, QuerySpan span,
    std::optional<std::uint64_t> diagonal_band_width)
    : ToDenseMatrix(std::make_shared<const PixelSelector>(std::move(sel)), n, span,
                    diagonal_band_width) {}

template <typename N, typename PixelSelector>
inline ToDenseMatrix<N, PixelSelector>::ToDenseMatrix(
    std::shared_ptr<const PixelSelector> sel, [[maybe_unused]] N n, QuerySpan span,
    std::optional<std::uint64_t> diagonal_band_width)
    : _sel(std::move(sel)), _span(span), _diagonal_band_width(diagonal_band_width) {
  if (!_sel) {
    throw std::runtime_error("hictk::transformers::ToDenseMatrix(): sel cannot be null");
  }
  if (chrom1() != chrom2() && span == QuerySpan::lower_triangle) {
    throw std::runtime_error(
        "hictk::transformers::ToDenseMatrix(): invalid parameters. Trans queries do not support "
        "span=QuerySpan::lower_triangle.");
  }

  validate_dtype();
}

template <typename N, typename PixelSelector>
inline auto ToDenseMatrix<N, PixelSelector>::operator()() -> MatrixT {
  if (!_sel) {
    return {};
  }

  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  const auto populate_upper_triangle =
      _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  auto matrix_setter = [](MatrixT& matrix, std::int64_t i1, std::int64_t i2, N count) noexcept {
    assert(i1 >= 0);
    assert(i1 < matrix.rows());
    assert(i2 >= 0);
    assert(i2 < matrix.cols());

    matrix(i1, i2) = count;
  };

  auto matrix = init_matrix();
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (chrom1() == chrom2() && _sel->coord1() != _sel->coord2()) {
      auto coord3 = _sel->coord1();
      auto coord4 = _sel->coord2();
      coord3.bin1 = std::min(coord3.bin1, coord4.bin1);
      coord3.bin2 = std::max(coord3.bin2, coord4.bin2);
      coord4 = coord3;

      const auto new_sel = _sel->fetch(coord3, coord4);
      internal::fill_matrix(new_sel.template begin<N>(), new_sel.template end<N>(),
                            _diagonal_band_width, internal::selector_is_symmetric_upper(new_sel),
                            matrix, matrix, matrix.rows(), matrix.cols(), row_offset(),
                            col_offset(), populate_lower_triangle, populate_upper_triangle,
                            matrix_setter);
      return matrix;
    }
  }

  internal::fill_matrix(_sel->template begin<N>(), _sel->template end<N>(), _diagonal_band_width,
                        internal::selector_is_symmetric_upper(*_sel), matrix, matrix, matrix.rows(),
                        matrix.cols(), row_offset(), col_offset(), populate_lower_triangle,
                        populate_upper_triangle, matrix_setter);
  return matrix;
}

template <typename N, typename PixelSelector>
inline std::string_view ToDenseMatrix<N, PixelSelector>::chrom1() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel->coord1().bin1.chrom().name();
  }
  return "all";
}

template <typename N, typename PixelSelector>
inline std::string_view ToDenseMatrix<N, PixelSelector>::chrom2() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel->coord2().bin1.chrom().name();
  }
  return "all";
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::num_bins(const PixelCoordinates& coords,
                                                              const BinTable& bins) noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    const auto span = coords.bin2.end() - coords.bin1.start();
    if (span == 0) {
      return static_cast<std::int64_t>(bins.size());
    }

    return static_cast<std::int64_t>(coords.bin2.id() - coords.bin1.id() + 1);
  }

  return static_cast<std::int64_t>(bins.size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::num_rows() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return num_bins(_sel->coord1(), _sel->bins());
  }
  return static_cast<std::int64_t>(_sel->bins().size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::num_cols() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return num_bins(_sel->coord2(), _sel->bins());
  }
  return static_cast<std::int64_t>(_sel->bins().size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::row_offset() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return offset(_sel->coord1());
  }
  return 0;
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::col_offset() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return offset(_sel->coord2());
  }
  return 0;
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::offset(
    const PixelCoordinates& coords) noexcept {
  constexpr auto bad_bin_id = std::numeric_limits<std::uint64_t>::max();
  return static_cast<std::int64_t>(coords.bin1.id() == bad_bin_id ? 0 : coords.bin1.id());
}

template <typename N, typename PixelSelector>
inline auto ToDenseMatrix<N, PixelSelector>::init_matrix() const -> MatrixT {
  assert(!!_sel);
  const auto& [weights1, weights2] = slice_weights(*_sel);

  if (weights1.size() == 0) {
    assert(weights2.size() == 0);
    return MatrixT::Zero(num_rows(), num_cols());
  }

  assert(weights2.size() != 0);
  return (weights1 * 0) * weights2.transpose();
}

template <typename N, typename PixelSelector>
inline std::pair<Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>,
                 Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>>
ToDenseMatrix<N, PixelSelector>::slice_weights(const cooler::PixelSelector& sel) const {
  return slice_weights(sel.weights(), sel.weights(), row_offset(), col_offset(), num_rows(),
                       num_cols());
}

template <typename N, typename PixelSelector>
inline std::pair<Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>,
                 Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>>
ToDenseMatrix<N, PixelSelector>::slice_weights(const hic::PixelSelector& sel) const {
  return slice_weights(
      sel.weights1(), sel.weights2(), static_cast<std::int64_t>(_sel->coord1().bin1.rel_id()),
      static_cast<std::int64_t>(_sel->coord2().bin1.rel_id()), num_rows(), num_cols());
}

template <typename N, typename PixelSelector>
inline std::pair<Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>,
                 Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>>
ToDenseMatrix<N, PixelSelector>::slice_weights(const hic::PixelSelectorAll& sel) const {
  return slice_weights(sel.weights(), sel.weights(), row_offset(), col_offset(), num_rows(),
                       num_cols());
}

template <typename N, typename PixelSelector>
inline std::pair<Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>,
                 Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>>
ToDenseMatrix<N, PixelSelector>::slice_weights(const hictk::PixelSelector& sel) const {
  return std::visit([&](const auto& sel_) { return slice_weights(sel_); }, sel.get());
}

template <typename N, typename PixelSelector>
inline std::pair<Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>,
                 Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor>>
ToDenseMatrix<N, PixelSelector>::slice_weights(const balancing::Weights& weights1,
                                               const balancing::Weights& weights2,
                                               std::int64_t offset1, std::int64_t offset2,
                                               std::int64_t size1, std::int64_t size2) {
  if constexpr (std::is_integral_v<N>) {
    return {};
  }

  if (weights1.empty() || weights2.empty()) {
    return {};
  }

  assert(offset1 + size1 <= static_cast<std::int64_t>(weights1.size()));
  assert(offset2 + size2 <= static_cast<std::int64_t>(weights2.size()));

  Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor> slice1(size1);

  for (std::int64_t i = 0; i < size1; ++i) {
    slice1(i) = conditional_static_cast<N>(weights1.at(static_cast<std::size_t>(offset1 + i),
                                                       balancing::Weights::Type::MULTIPLICATIVE));
  }

  const auto symmetric_query = &weights1 == &weights2 && offset1 == offset2 && size1 == size2;
  if (symmetric_query) {
    return std::make_pair(slice1, slice1);
  }

  Eigen::Matrix<N, Eigen::Dynamic, Eigen::RowMajor> slice2(size2);
  for (std::int64_t i = 0; i < size2; ++i) {
    slice2(i) = conditional_static_cast<N>(weights2.at(static_cast<std::size_t>(offset2 + i),
                                                       balancing::Weights::Type::MULTIPLICATIVE));
  }

  return std::make_pair(std::move(slice1), std::move(slice2));
}

template <typename N, typename PixelSelector>
inline void ToDenseMatrix<N, PixelSelector>::validate_dtype() const {
  assert(!!_sel);
  if constexpr (std::is_floating_point_v<N>) {
    return;
  }

  constexpr auto* msg =
      "hictk::transformers::ToDenseMatrix(): invalid parameters. n should be of floating-point "
      "type when fetching normalized interactions.";

  if constexpr (internal::has_weights_member_fx<PixelSelector>) {
    if (!_sel->weights().is_vector_of_ones()) {
      throw std::runtime_error(msg);
    }
  } else {
    if (!_sel->weights1().is_vector_of_ones() || !_sel->weights2().is_vector_of_ones()) {
      throw std::runtime_error(msg);
    }
  }
}

}  // namespace hictk::transformers
