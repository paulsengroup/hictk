// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cstdint>
#include <utility>

#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
inline ToDenseMatrix<N, PixelSelector>::ToDenseMatrix(PixelSelector&& sel, [[maybe_unused]] N n,
                                                      bool mirror)
    : _sel(std::move(sel)), _mirror(mirror) {}

template <typename N, typename PixelSelector>
inline auto ToDenseMatrix<N, PixelSelector>::operator()()
    -> Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic> {
  const auto offset1 = row_offset();
  const auto offset2 = col_offset();

  const auto num_rows_ = num_rows();
  const auto num_cols_ = num_cols();

  const auto mirror_matrix = _mirror && chrom1() == chrom2();

  using MatrixT = Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  MatrixT matrix = MatrixT::Zero(num_rows_, num_cols_);
  std::for_each(_sel.template begin<N>(), _sel.template end<N>(), [&](const ThinPixel<N>& p) {
    const auto i1 = static_cast<std::int64_t>(p.bin1_id - offset1);
    const auto i2 = static_cast<std::int64_t>(p.bin2_id - offset2);
    matrix(i1, i2) = p.count;

    if (mirror_matrix) {
      const auto delta = i2 - i1;
      if (delta >= 0 && delta < num_rows_ && i1 < num_cols_ && i2 < num_rows_) {
        matrix(i2, i1) = p.count;
      } else if ((delta < 0 || delta > num_cols_) && i1 < num_cols_ && i2 < num_rows_) {
        const auto i3 = static_cast<std::int64_t>(p.bin2_id - offset1);
        const auto i4 = static_cast<std::int64_t>(p.bin1_id - offset2);

        if (i3 >= 0 && i3 < num_rows_ && i4 >= 0 && i4 < num_cols_) {
          matrix(i3, i4) = p.count;
        }
      }
    }
  });
  return matrix;
}

template <typename N, typename PixelSelector>
inline std::int64_t ToDenseMatrix<N, PixelSelector>::num_rows() const noexcept {
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
inline std::int64_t ToDenseMatrix<N, PixelSelector>::num_cols() const noexcept {
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
inline std::uint64_t ToDenseMatrix<N, PixelSelector>::row_offset() const noexcept {
  constexpr auto bad_bin_id = std::numeric_limits<std::uint64_t>::max();
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel.coord1().bin1.id() == bad_bin_id ? 0 : _sel.coord1().bin1.id();
  }
  return 0;
}

template <typename N, typename PixelSelector>
inline std::uint64_t ToDenseMatrix<N, PixelSelector>::col_offset() const noexcept {
  constexpr auto bad_bin_id = std::numeric_limits<std::uint64_t>::max();
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    return _sel.coord2().bin1.id() == bad_bin_id ? 0 : _sel.coord2().bin1.id();
  }
  return 0;
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

}  // namespace hictk::transformers
