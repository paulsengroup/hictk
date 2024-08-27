// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Eigen/SparseCore>
#include <algorithm>
#include <cstdint>
#include <utility>

#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
inline ToSparseMatrix<N, PixelSelector>::ToSparseMatrix(PixelSelector&& sel, [[maybe_unused]] N n,
                                                        bool transpose)
    : _sel(std::move(sel)), _transpose(transpose) {}

template <typename N, typename PixelSelector>
inline auto ToSparseMatrix<N, PixelSelector>::operator()() -> Eigen::SparseMatrix<N> {
  const auto offset1 = row_offset();
  const auto offset2 = col_offset();

  auto num_rows_ = num_rows();
  auto num_cols_ = num_cols();

  if (_transpose) {
    std::swap(num_rows_, num_cols_);
  }

  Eigen::SparseMatrix<N> matrix(num_rows_, num_cols_);
  std::for_each(_sel.template begin<N>(), _sel.template end<N>(), [&](const ThinPixel<N>& p) {
    auto bin1 = static_cast<std::int64_t>(p.bin1_id - offset1);
    auto bin2 = static_cast<std::int64_t>(p.bin2_id - offset2);
    if (_transpose) {
      std::swap(bin1, bin2);
    }

    matrix.insert(bin1, bin2) = p.count;
  });
  matrix.makeCompressed();
  return matrix;
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::num_rows() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    const auto span = _sel.coord1().bin2.end() - _sel.coord1().bin1.start();
    if (span == 0) {
      return static_cast<std::int64_t>(_sel.bins().size());
    }

    const auto bin_size = _sel.bins().resolution();
    return static_cast<std::int64_t>((span + bin_size - 1) / bin_size);
  }

  return static_cast<std::int64_t>(_sel.bins().size());
}

template <typename N, typename PixelSelector>
inline std::int64_t ToSparseMatrix<N, PixelSelector>::num_cols() const noexcept {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    const auto span = _sel.coord2().bin2.end() - _sel.coord2().bin1.start();
    if (span == 0) {
      return static_cast<std::int64_t>(_sel.bins().size());
    }

    const auto bin_size = _sel.bins().resolution();
    return static_cast<std::int64_t>((span + bin_size - 1) / bin_size);
  }

  return static_cast<std::int64_t>(_sel.bins().size());
}

template <typename N, typename PixelSelector>
inline std::uint64_t ToSparseMatrix<N, PixelSelector>::row_offset() const noexcept {
  constexpr auto bad_bin_id = std::numeric_limits<std::uint64_t>::max();
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel.coord1().bin1.id() == bad_bin_id ? 0 : _sel.coord1().bin1.id();
  }
  return 0;
}

template <typename N, typename PixelSelector>
inline std::uint64_t ToSparseMatrix<N, PixelSelector>::col_offset() const noexcept {
  constexpr auto bad_bin_id = std::numeric_limits<std::uint64_t>::max();
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel.coord2().bin1.id() == bad_bin_id ? 0 : _sel.coord2().bin1.id();
  }
  return 0;
}

}  // namespace hictk::transformers
