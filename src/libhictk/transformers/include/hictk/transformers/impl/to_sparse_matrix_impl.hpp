// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Eigen/SparseCore>
#include <algorithm>
#include <cstdint>
#include <type_traits>
#include <utility>

#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
inline ToSparseMatrix<N, PixelSelector>::ToSparseMatrix(PixelSelector&& sel, [[maybe_unused]] N n)
    : _sel(std::move(sel)) {}

namespace internal {
template <typename T, typename = std::void_t<>>
inline constexpr bool has_coord1_member_fx = false;

template <typename T>
inline constexpr bool has_coord1_member_fx<T, std::void_t<decltype(std::declval<T>().coord1())>> =
    true;
}  // namespace internal

template <typename N, typename PixelSelector>
inline auto ToSparseMatrix<N, PixelSelector>::operator()() -> Eigen::SparseMatrix<N> {
  const auto offset1 = row_offset();
  const auto offset2 = col_offset();

  Eigen::SparseMatrix<N> matrix(num_rows(), num_cols());
  std::for_each(_sel.template begin<N>(), _sel.template end<N>(), [&](const ThinPixel<N>& p) {
    matrix.insert(static_cast<std::int64_t>(p.bin1_id - offset1),
                  static_cast<std::int64_t>(p.bin2_id - offset2)) = p.count;
  });
  matrix.makeCompressed();
  return matrix;
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::num_rows() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    const auto bin_size = _sel.bins().resolution();
    const auto span = _sel.coord1().bin2.end() - _sel.coord1().bin1.start();
    return static_cast<std::int64_t>((span + bin_size - 1) / bin_size);
  }

  return static_cast<std::int64_t>(_sel.bins().size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::num_cols() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    const auto bin_size = _sel.bins().resolution();
    const auto span = _sel.coord2().bin2.end() - _sel.coord2().bin1.start();
    return static_cast<std::int64_t>((span + bin_size - 1) / bin_size);
  }

  return static_cast<std::int64_t>(_sel.bins().size());
}

template <typename N, typename PixelSelector>
inline std::uint64_t ToSparseMatrix<N, PixelSelector>::row_offset() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel.coord1().bin1.id();
  }
  return 0;
}

template <typename N, typename PixelSelector>
inline std::uint64_t ToSparseMatrix<N, PixelSelector>::col_offset() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel.coord2().bin1.id();
  }
  return 0;
}

}  // namespace hictk::transformers
