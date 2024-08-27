// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cstdint>
#include <utility>

#include "hictk/pixel.hpp"
#include "hictk/transformers/impl/common.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
inline ToDenseMatrix<N, PixelSelector>::ToDenseMatrix(PixelSelector&& sel, [[maybe_unused]] N n,
                                                      bool mirror)
    : _sel(std::move(sel)), _mirror(mirror) {}

template <typename N, typename PixelSelector>
inline auto ToDenseMatrix<N, PixelSelector>::operator()() -> MatrixT {
  MatrixT matrix = MatrixT::Zero(num_rows(), num_cols());
  const auto mirror_matrix = _mirror && chrom1() == chrom2();

  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (chrom1() == chrom2() && _sel.coord2().bin2 > _sel.coord1().bin2) {
      auto coord3 = _sel.coord1();
      auto coord4 = _sel.coord2();
      coord3.bin1 = std::min(coord3.bin1, coord4.bin1);  // probably unnecessary
      coord3.bin2 = std::max(coord3.bin2, coord4.bin2);
      coord4 = coord3;

      fill_matrix(_sel.fetch(coord3, coord4), matrix, row_offset(), col_offset(), mirror_matrix);
      return matrix;
    }
  }

  fill_matrix(_sel, matrix, row_offset(), col_offset(), mirror_matrix);
  return matrix;
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
inline void ToDenseMatrix<N, PixelSelector>::fill_matrix(const PixelSelector& sel, MatrixT& buffer,
                                                         std::int64_t offset1, std::int64_t offset2,
                                                         bool mirror_matrix) {
  std::for_each(sel.template begin<N>(), sel.template end<N>(), [&](const ThinPixel<N>& p) {
    const auto i1 = static_cast<std::int64_t>(p.bin1_id) - offset1;
    const auto i2 = static_cast<std::int64_t>(p.bin2_id) - offset2;

    if (i1 >= 0 && i1 < buffer.rows() && i2 >= 0 && i2 < buffer.cols()) {
      buffer(i1, i2) = p.count;
    }

    if (mirror_matrix) {
      const auto i3 = static_cast<std::int64_t>(p.bin2_id) - offset1;
      const auto i4 = static_cast<std::int64_t>(p.bin1_id) - offset2;

      if (i3 >= 0 && i3 < buffer.rows() && i4 >= 0 && i4 < buffer.cols()) {
        buffer(i3, i4) = p.count;
      }
    }
  });
}

}  // namespace hictk::transformers
