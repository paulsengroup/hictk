// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Eigen/SparseCore>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <utility>

#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
inline ToSparseMatrix<N, PixelSelector>::ToSparseMatrix(PixelSelector sel, [[maybe_unused]] N n,
                                                        QuerySpan span)
    : ToSparseMatrix(std::make_shared<const PixelSelector>(std::move(sel)), n, span) {}

template <typename N, typename PixelSelector>
inline ToSparseMatrix<N, PixelSelector>::ToSparseMatrix(std::shared_ptr<const PixelSelector> sel,
                                                        [[maybe_unused]] N n, QuerySpan span)
    : _sel(std::move(sel)), _span(span) {
  if (!_sel) {
    throw std::runtime_error("hictk::transformers::ToSparseMatrix(): sel cannot be null");
  }
  if (chrom1() != chrom2() && span == QuerySpan::lower_triangle) {
    throw std::runtime_error(
        "hictk::transformers::ToSparseMatrix(): invalid parameters. Trans queries do not support "
        "span=QuerySpan::lower_triangle.");
  }
  validate_dtype();
}

template <typename N, typename PixelSelector>
inline auto ToSparseMatrix<N, PixelSelector>::operator()() -> MatrixT {
  assert(!!_sel);
  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  const auto populate_upper_triangle =
      _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  auto matrix_setter = [](MatrixT& matrix, std::int64_t i1, std::int64_t i2, N count) {
    matrix.insert(i1, i2) = count;
  };

  MatrixT matrix(num_rows(), num_cols());
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (chrom1() == chrom2() && _sel->coord1() != _sel->coord2()) {
      auto coord3 = _sel->coord1();
      auto coord4 = _sel->coord2();
      coord3.bin1 = std::min(coord3.bin1, coord4.bin1);
      coord3.bin2 = std::max(coord3.bin2, coord4.bin2);
      coord4 = coord3;

      const auto new_sel = _sel->fetch(coord3, coord4);
      internal::fill_matrix(new_sel.template begin<N>(), new_sel.template end<N>(),
                            internal::selector_is_symmetric_upper(new_sel), matrix, matrix.rows(),
                            matrix.cols(), row_offset(), col_offset(), populate_lower_triangle,
                            populate_upper_triangle, matrix_setter);
      matrix.makeCompressed();
      return matrix;
    }
  }

  internal::fill_matrix(_sel->template begin<N>(), _sel->template end<N>(),
                        internal::selector_is_symmetric_upper(*_sel), matrix, matrix.rows(),
                        matrix.cols(), row_offset(), col_offset(), populate_lower_triangle,
                        populate_upper_triangle, matrix_setter);
  matrix.makeCompressed();
  return matrix;
}

template <typename N, typename PixelSelector>
inline std::string_view ToSparseMatrix<N, PixelSelector>::chrom1() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel->coord1().bin1.chrom().name();
  }
  return "all";
}

template <typename N, typename PixelSelector>
inline std::string_view ToSparseMatrix<N, PixelSelector>::chrom2() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel->coord2().bin1.chrom().name();
  }
  return "all";
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::num_bins(const PixelCoordinates& coords,
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
inline std::int64_t ToSparseMatrix<N, PixelSelector>::num_rows() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return num_bins(_sel->coord1(), _sel->bins());
  }
  return static_cast<std::int64_t>(_sel->bins().size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::num_cols() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return num_bins(_sel->coord2(), _sel->bins());
  }
  return static_cast<std::int64_t>(_sel->bins().size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::row_offset() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return offset(_sel->coord1());
  }
  return 0;
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::col_offset() const noexcept {
  assert(!!_sel);
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return offset(_sel->coord2());
  }
  return 0;
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::offset(
    const PixelCoordinates& coords) noexcept {
  constexpr auto bad_bin_id = std::numeric_limits<std::uint64_t>::max();
  return static_cast<std::int64_t>(coords.bin1.id() == bad_bin_id ? 0 : coords.bin1.id());
}

template <typename N, typename PixelSelector>
inline void ToSparseMatrix<N, PixelSelector>::validate_dtype() const {
  assert(!!_sel);
  if constexpr (std::is_floating_point_v<N>) {
    return;
  }

  constexpr auto* msg =
      "hictk::transformers::ToSparseMatrix(): invalid parameters. n should be of floating-point "
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
