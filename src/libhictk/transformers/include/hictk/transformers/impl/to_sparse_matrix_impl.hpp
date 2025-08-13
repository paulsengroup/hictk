// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#if __has_include(<eigen3/Eigen/SparseCore>)
#include <eigen3/Eigen/SparseCore>
#else
#include <Eigen/SparseCore>
#endif

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/pixel.hpp"
#include "hictk/transformers/common.hpp"
#include "hictk/transformers/diagonal_band.hpp"

namespace hictk::transformers {

template <typename N, typename PixelSelector>
inline ToSparseMatrix<N, PixelSelector>::ToSparseMatrix(
    PixelSelector sel, [[maybe_unused]] N n, QuerySpan span, bool minimize_memory_usage,
    std::optional<std::uint64_t> diagonal_band_width)
    : ToSparseMatrix(std::make_shared<const PixelSelector>(std::move(sel)), n, span,
                     minimize_memory_usage, diagonal_band_width) {}

template <typename N, typename PixelSelector>
inline ToSparseMatrix<N, PixelSelector>::ToSparseMatrix(
    std::shared_ptr<const PixelSelector> sel, [[maybe_unused]] N n, QuerySpan span,
    bool minimize_memory_usage, std::optional<std::uint64_t> diagonal_band_width)
    : _sel(std::move(sel)),
      _span(span),
      _minimize_memory_usage(minimize_memory_usage),
      _diagonal_band_width(diagonal_band_width) {
  if (!_sel) {
    throw std::runtime_error("hictk::transformers::ToSparseMatrix(): sel cannot be null");
  }
  if (chrom1() != chrom2() && _span == QuerySpan::lower_triangle) {
    throw std::runtime_error(
        "hictk::transformers::ToSparseMatrix(): invalid parameters. Trans queries do not support "
        "span=QuerySpan::lower_triangle.");
  }
  validate_dtype();
}

template <typename N, typename PixelSelector>
inline auto ToSparseMatrix<N, PixelSelector>::operator()() -> MatrixT {
  if (!_sel) {
    return {};
  }

  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  const auto populate_upper_triangle =
      _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (chrom1() == chrom2() && _sel->coord1() != _sel->coord2()) {
      auto coord3 = _sel->coord1();
      auto coord4 = _sel->coord2();
      coord3.bin1 = std::min(coord3.bin1, coord4.bin1);
      coord3.bin2 = std::max(coord3.bin2, coord4.bin2);
      coord4 = coord3;

      if (_minimize_memory_usage) {
        return fill_matrix_low_mem(_sel->fetch(coord3, coord4), populate_upper_triangle,
                                   populate_lower_triangle);
      }
      return fill_matrix_fast(_sel->fetch(coord3, coord4), populate_upper_triangle,
                              populate_lower_triangle);
    }
  }

  if (_minimize_memory_usage) {
    return fill_matrix_low_mem(*_sel, populate_upper_triangle, populate_lower_triangle);
  }
  return fill_matrix_fast(*_sel, populate_upper_triangle, populate_lower_triangle);
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

template <typename N, typename PixelSelector>
template <typename PixelIt_>
inline auto ToSparseMatrix<N, PixelSelector>::fill_matrix_fast(PixelIt_ first, const PixelIt_& last,
                                                               bool selector_is_symmetric_upper,
                                                               bool populate_upper_triangle,
                                                               bool populate_lower_triangle) const
    -> MatrixT {
  MatrixRowMajor matrix_ut(num_rows(), num_cols());
  MatrixColMajor matrix_lt{};

  const auto mirror_interactions = interactions_can_be_mirrored();
  auto transpose_interactions = interactions_can_be_transposed();

  if (_span == QuerySpan::full && transpose_interactions) {
    // output matrix will be computed as M = MU + MU.T, where MU is row-major
    populate_upper_triangle = true;
    populate_lower_triangle = false;
  } else if (selector_is_symmetric_upper && mirror_interactions) {
    // output matrix will be computed as M = MU + ML, where MU is row-major and ML is col-major
    matrix_lt = MatrixColMajor(num_rows(), num_cols());
    transpose_interactions = false;
  }

  std::vector<ThinPixel<N>> row_buff{};
  std::int64_t reserved_size_ut{0};
  std::int64_t reserved_size_lt{0};

  if (matrix_lt.size() == 0) {
    while (first != last) {
      fill_row(first, last, matrix_ut, matrix_ut, reserved_size_ut, reserved_size_lt, row_buff,
               selector_is_symmetric_upper, row_offset(), col_offset(), populate_lower_triangle,
               populate_upper_triangle);
    }
  } else {
    while (first != last) {
      fill_row(first, last, matrix_ut, matrix_lt, reserved_size_lt, reserved_size_lt, row_buff,
               selector_is_symmetric_upper, row_offset(), col_offset(), populate_lower_triangle,
               populate_upper_triangle);
    }
  }

  if (matrix_lt.nonZeros() != 0) {
    matrix_ut += MatrixRowMajor(std::move(matrix_lt));
  } else if (_span == QuerySpan::full && transpose_interactions) {
    matrix_ut +=
        MatrixRowMajor(matrix_ut.template triangularView<Eigen::StrictlyUpper>().transpose());
  }

  matrix_ut.makeCompressed();

  [[maybe_unused]] const auto space_overhead =
      1.0 - (static_cast<double>(matrix_ut.nonZeros()) /
             static_cast<double>(reserved_size_ut + reserved_size_lt));

  SPDLOG_DEBUG(FMT_STRING("ToSparseMatrix::fill_matrix_fast(): space overhead: {:.3f}% ({} nnz)"),
               100.0 * space_overhead,
               (reserved_size_ut + reserved_size_lt) - matrix_ut.nonZeros());

  return matrix_ut;
}

template <typename N, typename PixelSelector>
inline auto ToSparseMatrix<N, PixelSelector>::fill_matrix_fast(const PixelSelector& sel,
                                                               bool populate_upper_triangle,
                                                               bool populate_lower_triangle) const
    -> MatrixT {
  const auto selector_is_symmetric_upper = internal::selector_is_symmetric_upper(sel);
  if (HICTK_UNLIKELY(_diagonal_band_width.has_value())) {
    const DiagonalBand band_sel(sel.template begin<N>(), sel.template end<N>(),
                                *_diagonal_band_width);
    return fill_matrix_fast(band_sel.begin(), band_sel.end(), selector_is_symmetric_upper,
                            populate_upper_triangle, populate_lower_triangle);
  }

  return fill_matrix_fast(sel.template begin<N>(), sel.template end<N>(),
                          selector_is_symmetric_upper, populate_upper_triangle,
                          populate_lower_triangle);
}
template <typename N, typename PixelSelector>
inline auto ToSparseMatrix<N, PixelSelector>::fill_matrix_low_mem(
    const PixelSelector& sel, bool populate_upper_triangle, bool populate_lower_triangle) const
    -> MatrixT {
  auto matrix_setter = [](auto& m, std::int64_t i1, std::int64_t i2, N count) noexcept {
    assert(i1 >= 0);
    assert(i1 < m.rows());
    assert(i2 >= 0);
    assert(i2 < m.cols());

    m.insert(i1, i2) = count;
  };

  const auto selector_is_symmetric_upper = internal::selector_is_symmetric_upper(sel);
  auto matrix = pre_allocate_matrix(sel, populate_upper_triangle, populate_lower_triangle);

  internal::fill_matrix(sel.template begin<N>(), sel.template end<N>(), _diagonal_band_width,
                        selector_is_symmetric_upper, matrix, matrix, matrix.rows(), matrix.cols(),
                        row_offset(), col_offset(), populate_lower_triangle,
                        populate_upper_triangle, matrix_setter);

  matrix.makeCompressed();

  return matrix;
}

template <typename N, typename PixelSelector>
inline bool ToSparseMatrix<N, PixelSelector>::interactions_can_be_transposed() const noexcept {
  if (!_sel) {
    return false;
  }

  const auto selector_is_symmetric_upper = internal::selector_is_symmetric_upper(*_sel);
  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;

  // NOLINTNEXTLINE(misc-const-correctness)
  bool res = selector_is_symmetric_upper && populate_lower_triangle &&
             !internal::has_coord1_member_fx<PixelSelector>;
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (selector_is_symmetric_upper) {
      res = populate_lower_triangle && _sel->coord1() == _sel->coord2();
    }
  }

  return res;
}

template <typename N, typename PixelSelector>
inline bool ToSparseMatrix<N, PixelSelector>::interactions_can_be_mirrored() const noexcept {
  if (!_sel) {
    return false;
  }

  const auto selector_is_symmetric_upper = internal::selector_is_symmetric_upper(*_sel);
  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;

  // NOLINTNEXTLINE(misc-const-correctness)
  bool res = selector_is_symmetric_upper && populate_lower_triangle &&
             !internal::has_coord1_member_fx<PixelSelector>;
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (selector_is_symmetric_upper) {
      res = populate_lower_triangle && _sel->coord1().bin1.chrom() == _sel->coord2().bin1.chrom();
    }
  }

  return res;
}

template <typename N, typename PixelSelector>
template <typename Matrix, typename PixelIt_>
inline void ToSparseMatrix<N, PixelSelector>::fill_row(
    PixelIt_& first_pixel, const PixelIt_& last_pixel, MatrixRowMajor& matrix_ut, Matrix& matrix_lt,
    std::int64_t& reserved_size_ut, std::int64_t& reserved_size_lt,
    std::vector<ThinPixel<N>>& buffer, bool symmetric_upper, std::int64_t offset1,
    std::int64_t offset2, bool populate_lower_triangle, bool populate_upper_triangle) {
  buffer.clear();
  if (first_pixel == last_pixel) {
    return;
  }

  std::int64_t num_pixels_ut{};

  const auto row = first_pixel->bin1_id;
  while (first_pixel != last_pixel) {
    auto pixel = *first_pixel;
    if (pixel.bin1_id != row) {
      break;
    }

    num_pixels_ut +=
        !symmetric_upper || (populate_upper_triangle && pixel.bin2_id >= pixel.bin1_id);
    buffer.emplace_back(std::move(pixel));
    std::ignore = ++first_pixel;
  }

  const auto size_ut_upper_bound = matrix_ut.rows() * num_pixels_ut;
  if (size_ut_upper_bound > reserved_size_ut) {
    const auto new_size =
        static_cast<std::int64_t>(1.25 * static_cast<double>(size_ut_upper_bound));
    SPDLOG_DEBUG(FMT_STRING("ToSparseMatrix::fill_row(): resizing UT matrix from {} to {}..."),
                 reserved_size_ut, new_size);
    matrix_ut.reserve(std::vector(static_cast<std::size_t>(matrix_ut.rows()),
                                  (new_size + matrix_ut.rows() - 1) / matrix_ut.rows()));
    reserved_size_ut = new_size;
  }

  if constexpr (!std::is_same_v<MatrixRowMajor, Matrix>) {
    const auto num_pixels_lt = static_cast<std::int64_t>(buffer.size()) - num_pixels_ut;
    const auto size_lt_upper_bound = matrix_lt.cols() * num_pixels_lt;
    if (size_lt_upper_bound > reserved_size_lt) {
      const auto new_size =
          static_cast<std::int64_t>(1.25 * static_cast<double>(size_lt_upper_bound));
      SPDLOG_DEBUG(FMT_STRING("ToSparseMatrix::fill_row(): resizing LT matrix from {} to {}..."),
                   reserved_size_lt, new_size);
      matrix_lt.reserve(std::vector(static_cast<std::size_t>(matrix_lt.cols()),
                                    (new_size + matrix_lt.cols() - 1) / matrix_lt.cols()));
      reserved_size_lt = new_size;
    }
  }

  auto matrix_setter = [](auto& m, std::int64_t i1, std::int64_t i2, N count) noexcept {
    assert(i1 >= 0);
    assert(i1 < m.rows());
    assert(i2 >= 0);
    assert(i2 < m.cols());

    m.insert(i1, i2) = count;
  };

  internal::fill_matrix(buffer.begin(), buffer.end(), std::nullopt, symmetric_upper, matrix_ut,
                        matrix_lt, matrix_ut.rows(), matrix_ut.cols(), offset1, offset2,
                        populate_lower_triangle, populate_upper_triangle, matrix_setter);
}

template <typename N, typename PixelSelector>
inline auto ToSparseMatrix<N, PixelSelector>::pre_allocate_matrix(
    const PixelSelector& sel, bool populate_upper_triangle, bool populate_lower_triangle) const
    -> MatrixT {
  auto setter = [](std::vector<std::int64_t>& buff, std::int64_t i1,
                   [[maybe_unused]] std::int64_t i2,
                   [[maybe_unused]] std::int_fast8_t count) noexcept {
    assert(i1 >= 0);
    assert(static_cast<std::size_t>(i1) < buff.size());
    ++buff[static_cast<std::size_t>(i1)];
  };

  const auto selector_is_symmetric_upper = internal::selector_is_symmetric_upper(sel);

  std::vector<std::int64_t> nnzs(static_cast<std::size_t>(num_rows()), 0);
  internal::fill_matrix(
      sel.template begin<std::int_fast8_t>(), sel.template end<std::int_fast8_t>(),
      _diagonal_band_width, selector_is_symmetric_upper, nnzs, nnzs, num_rows(), num_cols(),
      row_offset(), col_offset(), populate_lower_triangle, populate_upper_triangle, setter);

  MatrixT matrix(num_rows(), num_cols());
  matrix.reserve(nnzs);

  return matrix;
}

}  // namespace hictk::transformers
