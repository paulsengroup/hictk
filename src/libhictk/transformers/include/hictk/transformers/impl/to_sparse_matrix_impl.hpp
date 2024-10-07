// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#include <Eigen/SparseCore>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

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
  auto populate_lower_triangle = _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  auto populate_upper_triangle = _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  using SparseMatrix = std::variant<MatrixRowMajor, MatrixColMajor>;

  const auto selector_is_symmetric_upper = internal::selector_is_symmetric_upper(*_sel);
  auto matrix_ut = std::make_shared<SparseMatrix>(MatrixRowMajor(num_rows(), num_cols()));
  auto matrix_lt = matrix_ut;

  bool transpose_matrix = selector_is_symmetric_upper && !populate_upper_triangle &&
                          populate_lower_triangle && !internal::has_coord1_member_fx<PixelSelector>;
  bool mirror_matrix = selector_is_symmetric_upper && populate_upper_triangle &&
                       populate_lower_triangle && !internal::has_coord1_member_fx<PixelSelector>;
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (selector_is_symmetric_upper) {
      transpose_matrix = !populate_upper_triangle && populate_lower_triangle &&
                         chrom1() == chrom2() && _sel->coord1() == _sel->coord2();
      mirror_matrix = populate_upper_triangle && populate_lower_triangle && chrom1() == chrom2() &&
                      _sel->coord1() == _sel->coord2();

      if (transpose_matrix && mirror_matrix) {
        populate_upper_triangle = true;
        populate_lower_triangle = false;
      } else if (transpose_matrix) {
        populate_upper_triangle = true;
        populate_lower_triangle = false;
      } else if (mirror_matrix) {
        matrix_lt = std::make_shared<SparseMatrix>(MatrixColMajor(num_rows(), num_cols()));
      }
    }
  }

  std::vector<ThinPixel<N>> row_buff{};
  std::int64_t reserved_size_ut{0};
  std::int64_t reserved_size_lt{0};

  auto first_pixel = _sel->template begin<N>();
  auto last_pixel = _sel->template end<N>();
  std::unique_ptr<const PixelSelector> sel2{};
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (chrom1() == chrom2() && _sel->coord1() != _sel->coord2()) {
      auto coord3 = _sel->coord1();
      auto coord4 = _sel->coord2();
      coord3.bin1 = std::min(coord3.bin1, coord4.bin1);
      coord3.bin2 = std::max(coord3.bin2, coord4.bin2);
      coord4 = coord3;

      sel2 = std::make_unique<const PixelSelector>(_sel->fetch(coord3, coord4));
      first_pixel = sel2->template begin<N>();
      last_pixel = sel2->template end<N>();
    }
  }

  auto matrix_setter = [](auto& matrix, std::int64_t i1, std::int64_t i2, N count) {
    matrix.insert(i1, i2) = count;
  };

  auto& matrix1 = std::get<MatrixRowMajor>(*matrix_ut);
  std::visit(
      [&](auto& matrix2) {
        while (first_pixel != last_pixel) {
          first_pixel = fill_row(std::move(first_pixel), last_pixel, matrix1, matrix2,
                                 reserved_size_ut, reserved_size_lt, row_buff,
                                 selector_is_symmetric_upper, row_offset(), col_offset(),
                                 populate_lower_triangle, populate_upper_triangle, matrix_setter);
        }
      },
      *matrix_lt);

  if (mirror_matrix && transpose_matrix) {
    assert(matrix_ut == matrix_lt);
    matrix1 += MatrixRowMajor(matrix1.template triangularView<Eigen::StrictlyUpper>().transpose());
  } else if (transpose_matrix) {
    assert(matrix_ut == matrix_lt);
    matrix1 = MatrixRowMajor(matrix1.transpose());
  } else if (mirror_matrix) {
    assert(matrix_ut != matrix_lt);
    matrix1 += MatrixRowMajor(std::get<MatrixColMajor>(*matrix_lt));
  }

  matrix1.makeCompressed();

  [[maybe_unused]] const auto space_overhead =
      1.0 - (static_cast<double>(matrix1.nonZeros()) /
             static_cast<double>(reserved_size_ut + reserved_size_lt));

  SPDLOG_DEBUG(FMT_STRING("ToSparseMatrix::operator()(): space overhead: {:.3f}% ({} nnz)"),
               100.0 * space_overhead, (reserved_size_ut + reserved_size_lt) - matrix1.nonZeros());
  return matrix1;
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
template <typename Matrix, typename SetterOp>
inline auto ToSparseMatrix<N, PixelSelector>::fill_row(
    PixelIt first_pixel, PixelIt last_pixel, MatrixRowMajor& matrix_ut, Matrix& matrix_lt,
    std::int64_t& reserved_size_ut, std::int64_t& reserved_size_lt,
    std::vector<ThinPixel<N>>& buffer, bool symmetric_upper, std::int64_t offset1,
    std::int64_t offset2, bool populate_lower_triangle, bool populate_upper_triangle,
    SetterOp matrix_setter) -> PixelIt {
  buffer.clear();
  if (first_pixel == last_pixel) {
    return last_pixel;
  }

  std::int64_t num_pixels_ut{};

  const auto row = first_pixel->bin1_id;
  while (first_pixel != last_pixel) {
    auto pixel = *first_pixel;
    if (pixel.bin1_id != row) {
      break;
    }
    num_pixels_ut += !symmetric_upper || pixel.bin2_id >= pixel.bin1_id;
    buffer.emplace_back(std::move(pixel));
    std::ignore = ++first_pixel;
  }

  const auto size_ut_upper_bound = matrix_ut.rows() * num_pixels_ut;
  if (size_ut_upper_bound > reserved_size_ut) {
    const auto new_size =
        static_cast<std::int64_t>(1.25 * static_cast<double>(size_ut_upper_bound));
    SPDLOG_DEBUG(
        FMT_STRING("ToSparseMatrix::operator()(): resizing UT sparse matrix from {} to {}..."),
        reserved_size_ut, new_size);
    matrix_ut.reserve(std::vector(static_cast<std::size_t>(matrix_ut.rows()),
                                  (new_size + matrix_ut.rows() - 1) / matrix_ut.rows()));
    reserved_size_ut = new_size;
  }

  if constexpr (!std::is_same_v<MatrixRowMajor, Matrix>) {
    const auto num_pixels_lt = static_cast<std::int64_t>(buffer.size()) - num_pixels_ut;
    const auto size_lt_upper_bound = matrix_lt.rows() * num_pixels_lt;
    if (size_lt_upper_bound > reserved_size_lt) {
      const auto new_size =
          static_cast<std::int64_t>(1.25 * static_cast<double>(size_lt_upper_bound));
      SPDLOG_DEBUG(
          FMT_STRING("ToSparseMatrix::operator()(): resizing LT sparse matrix from {} to {}..."),
          reserved_size_lt, new_size);
      matrix_lt.reserve(std::vector(static_cast<std::size_t>(matrix_lt.rows()),
                                    (new_size + matrix_lt.rows() - 1) / matrix_lt.rows()));
      reserved_size_lt = new_size;
    }
  }

  internal::fill_matrix(buffer.begin(), buffer.end(), symmetric_upper, matrix_ut, matrix_lt,
                        matrix_ut.rows(), matrix_ut.cols(), offset1, offset2,
                        populate_lower_triangle, populate_upper_triangle, matrix_setter);
  return first_pixel;
}

}  // namespace hictk::transformers
