// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <utility>
#include <variant>

#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
inline ToDenseMatrix<N, PixelSelector>::ToDenseMatrix(PixelSelector&& sel, [[maybe_unused]] N n,
                                                      QuerySpan span)
    : _sel(std::move(sel)), _span(span) {
  if (chrom1() != chrom2() && span == QuerySpan::lower_triangle) {
    throw std::runtime_error(
        "hictk::transformers::ToDenseMatrix(): invalid parameters. Trans queries do not support "
        "span=QuerySpan::lower_triangle.");
  }
}

template <typename N, typename PixelSelector>
inline auto ToDenseMatrix<N, PixelSelector>::operator()() -> MatrixT {
  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  const auto populate_upper_triangle =
      _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  auto matrix_setter = [](MatrixT& matrix, std::int64_t i1, std::int64_t i2, N count) {
    matrix(i1, i2) = count;
  };

  MatrixT matrix = MatrixT::Zero(num_rows(), num_cols());
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (chrom1() == chrom2() && _sel.coord1() != _sel.coord2()) {
      auto coord3 = _sel.coord1();
      auto coord4 = _sel.coord2();
      coord3.bin1 = std::min(coord3.bin1, coord4.bin1);
      coord3.bin2 = std::max(coord3.bin2, coord4.bin2);
      coord4 = coord3;

      const auto new_sel = _sel.fetch(coord3, coord4);
      internal::fill_matrix<N>(new_sel, matrix, matrix.rows(), matrix.cols(), row_offset(),
                               col_offset(), populate_lower_triangle, populate_upper_triangle,
                               matrix_setter);
      mask_bad_bins(new_sel, matrix);
      return matrix;
    }
  }

  internal::fill_matrix<N>(_sel, matrix, matrix.rows(), matrix.cols(), row_offset(), col_offset(),
                           populate_lower_triangle, populate_upper_triangle, matrix_setter);
  mask_bad_bins(_sel, matrix);
  return matrix;
}

template <typename N, typename PixelSelector>
inline std::string_view ToDenseMatrix<N, PixelSelector>::chrom1() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel.coord1().bin1.chrom().name();
  }
  return "all";
}

template <typename N, typename PixelSelector>
inline std::string_view ToDenseMatrix<N, PixelSelector>::chrom2() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel.coord2().bin1.chrom().name();
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
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return num_bins(_sel.coord1(), _sel.bins());
  }
  return static_cast<std::int64_t>(_sel.bins().size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::num_cols() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return num_bins(_sel.coord2(), _sel.bins());
  }
  return static_cast<std::int64_t>(_sel.bins().size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::row_offset() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return offset(_sel.coord1());
  }
  return 0;
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::col_offset() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return offset(_sel.coord2());
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
inline void ToDenseMatrix<N, PixelSelector>::mask_bad_bins(const hictk::PixelSelector& sel,
                                                           MatrixT& buffer) {
  std::visit([&](const auto& sel_) { mask_bad_bins(sel_, buffer); }, sel.get());
}

template <typename N, typename PixelSelector>
inline void ToDenseMatrix<N, PixelSelector>::mask_bad_bins(const cooler::PixelSelector& sel,
                                                           MatrixT& buffer) {
  if constexpr (std::is_integral_v<N>) {
    return;
  }

  if (!sel.weights()) {
    return;
  }

  const auto offset1 = row_offset();
  const auto offset2 = col_offset();

  const auto weights = (*sel.weights())(balancing::Weights::Type::MULTIPLICATIVE);

  for (std::int64_t i = 0; i < buffer.rows(); ++i) {
    if (!std::isfinite(weights[static_cast<std::size_t>(offset1 + i)])) {
      buffer.row(i).setConstant(std::numeric_limits<N>::quiet_NaN());
    }
  }

  for (std::int64_t j = 0; j < buffer.cols(); ++j) {
    if (!std::isfinite(weights[static_cast<std::size_t>(offset2 + j)])) {
      buffer.col(j).setConstant(std::numeric_limits<N>::quiet_NaN());
    }
  }
}

template <typename N, typename PixelSelector>
inline void ToDenseMatrix<N, PixelSelector>::mask_bad_bins(const hic::PixelSelector& sel,
                                                           MatrixT& buffer) {
  if constexpr (std::is_integral_v<N>) {
    return;
  }

  if (sel.weights1()) {
    const auto chrom_offset = sel.bins().at(sel.coord1().bin1.chrom()).id();
    const auto offset = row_offset() - static_cast<std::int64_t>(chrom_offset);

    const auto& weights = sel.weights1();
    assert(weights.type() == balancing::Weights::Type::DIVISIVE);

    for (std::int64_t i = 0; i < buffer.rows(); ++i) {
      const auto w = weights[static_cast<std::size_t>(offset + i)];
      if (!std::isfinite(w) || w == 0) {
        buffer.row(i).setConstant(std::numeric_limits<N>::quiet_NaN());
      }
    }
  }

  if (sel.weights2()) {
    const auto chrom_offset = sel.bins().at(sel.coord2().bin1.chrom()).id();
    const auto offset = row_offset() - static_cast<std::int64_t>(chrom_offset);

    const auto& weights = sel.weights2();
    assert(weights.type() == balancing::Weights::Type::DIVISIVE);

    for (std::int64_t j = 0; j < buffer.cols(); ++j) {
      const auto w = weights[static_cast<std::size_t>(offset + j)];
      if (!std::isfinite(w) || w == 0) {
        buffer.col(j).setConstant(std::numeric_limits<N>::quiet_NaN());
      }
    }
  }
}

template <typename N, typename PixelSelector>
inline void ToDenseMatrix<N, PixelSelector>::mask_bad_bins(const hic::PixelSelectorAll& sel,
                                                           MatrixT& buffer) {
  if constexpr (std::is_integral_v<N>) {
    return;
  }

  const auto weights = sel.weights();
  assert(weights.type() == balancing::Weights::Type::DIVISIVE);

  for (std::int64_t i = 0; i < buffer.rows(); ++i) {
    const auto w = weights[static_cast<std::size_t>(i)];
    if (!std::isfinite(w) || w == 0) {
      buffer.row(i).setConstant(std::numeric_limits<N>::quiet_NaN());
    }
  }

  for (std::int64_t j = 0; j < buffer.cols(); ++j) {
    const auto w = weights[static_cast<std::size_t>(j)];
    if (!std::isfinite(w) || w == 0) {
      buffer.col(j).setConstant(std::numeric_limits<N>::quiet_NaN());
    }
  }
}

}  // namespace hictk::transformers
