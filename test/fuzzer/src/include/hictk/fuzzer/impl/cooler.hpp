// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/pybind11.h>

#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <string_view>

#include "hictk/cooler/cooler.hpp"
#include "hictk/fuzzer/common.hpp"
#include "hictk/reference.hpp"

namespace hictk::fuzzer::cooler {

template <typename N>
inline std::size_t COODataFrame<N>::size() const noexcept {
  assert(bin1_id.size() == bin2_id.size());
  assert(bin1_id.size() == count.size());

  return static_cast<std::size_t>(count.size());
}

template <typename N>
inline COODataFrame<N>& COODataFrame<N>::operator=(pybind11::object df) {
  bin1_id = pybind11::cast<NumpyArray<std::int64_t>>(df.attr("__getitem__")("bin1_id"));
  bin2_id = pybind11::cast<NumpyArray<std::int64_t>>(df.attr("__getitem__")("bin2_id"));
  if (pybind11::cast<bool>(df.attr("__contains__")("balanced"))) {
    if constexpr (std::is_integral_v<N>) {
      throw std::runtime_error(
          "fetching balanced interactions requires COODataFrame to be of floating-point type");
    }
    count = pybind11::cast<NumpyArray<N>>(df.attr("__getitem__")("balanced"));
  } else {
    count = pybind11::cast<NumpyArray<N>>(df.attr("__getitem__")("count"));
  }

  return *this;
}

template <typename N>
inline std::vector<ThinPixel<N>> COODataFrame<N>::to_vector() const {
  std::vector<ThinPixel<N>> buffer(size());
  to_vector(buffer);
  return buffer;
}

template <typename N>
void COODataFrame<N>::to_vector(std::vector<ThinPixel<N>>& buffer) const {
  buffer.resize(size());

  for (std::size_t i = 0; i < size(); ++i) {
    buffer[i] = ThinPixel<N>{static_cast<std::uint64_t>(bin1_id.at(i)),
                             static_cast<std::uint64_t>(bin2_id.at(i)), count.at(i)};
  }
}

template <typename N>
inline std::size_t BG2DataFrame<N>::size() const noexcept {
  assert(static_cast<std::int64_t>(chrom1.size()) == start1.size());
  assert(static_cast<std::int64_t>(chrom1.size()) == end1.size());
  assert(chrom1.size() == chrom2.size());
  assert(static_cast<std::int64_t>(chrom1.size()) == start2.size());
  assert(static_cast<std::int64_t>(chrom1.size()) == end2.size());
  assert(static_cast<std::int64_t>(chrom1.size()) == count.size());

  return static_cast<std::size_t>(count.size());
}

template <typename N>
inline BG2DataFrame<N>& BG2DataFrame<N>::operator=(pybind11::object df) {
  chrom1 = pybind11::cast<pybind11::list>(df.attr("__getitem__")("chrom1").attr("tolist")());
  start1 = pybind11::cast<NumpyArray<std::int64_t>>(df.attr("__getitem__")("start1"));
  end1 = pybind11::cast<NumpyArray<std::int64_t>>(df.attr("__getitem__")("end1"));
  chrom2 = pybind11::cast<pybind11::list>(df.attr("__getitem__")("chrom2").attr("tolist")());
  start2 = pybind11::cast<NumpyArray<std::int64_t>>(df.attr("__getitem__")("start2"));
  end2 = pybind11::cast<NumpyArray<std::int64_t>>(df.attr("__getitem__")("end2"));

  if (pybind11::cast<bool>(df.attr("__contains__")("balanced"))) {
    if constexpr (std::is_integral_v<N>) {
      throw std::runtime_error(
          "fetching balanced interactions requires BG2DataFrame to be of floating-point type");
    }
    count = pybind11::cast<NumpyArray<N>>(df.attr("__getitem__")("balanced"));
  } else {
    count = pybind11::cast<NumpyArray<N>>(df.attr("__getitem__")("count"));
  }

  return *this;
}

template <typename N>
inline std::vector<Pixel<N>> BG2DataFrame<N>::to_vector(const Reference& chroms) const {
  std::vector<Pixel<N>> buffer(size());
  to_vector(chroms, buffer);
  return buffer;
}

template <typename N>
void BG2DataFrame<N>::to_vector(const Reference& chroms, std::vector<Pixel<N>>& buffer) const {
  buffer.resize(size());

  auto chrom1_id = chrom1.begin();
  auto chrom2_id = chrom2.begin();
  for (std::size_t i = 0; i < size(); ++i) {
    buffer[i] = Pixel{chroms.at(pybind11::cast<std::string>(*chrom1_id)),
                      static_cast<std::uint32_t>(start1.at(i)),
                      static_cast<std::uint32_t>(end1.at(i)),
                      chroms.at(pybind11::cast<std::string>(*chrom2_id)),
                      static_cast<std::uint32_t>(start2.at(i)),
                      static_cast<std::uint32_t>(end2.at(i)),
                      count.at(i)};
    ++chrom1_id;
    ++chrom2_id;
  }
}

template <typename N>
inline EigenSparse<N> scipy_coo_to_eigen(pybind11::object obj) {
  const auto rows = pybind11::cast<NumpyArray<std::int64_t>>(obj.attr("row"));
  const auto cols = pybind11::cast<NumpyArray<std::int64_t>>(obj.attr("col"));
  const auto count = pybind11::cast<NumpyArray<N>>(obj.attr("data"));

  const auto num_rows = obj.attr("shape").attr("__getitem__")(0).cast<std::int64_t>();
  const auto num_cols = obj.attr("shape").attr("__getitem__")(1).cast<std::int64_t>();

  EigenSparse<N> m(num_rows, num_cols);
  if (rows.size() != 0) {
    const auto max_nnz_row = obj.attr("getnnz")(0).attr("max")().cast<std::int64_t>();
    m.reserve(Eigen::Matrix<std::int64_t, Eigen::Dynamic, 1>::Constant(num_cols, max_nnz_row));

    for (std::int64_t i = 0; i < rows.size(); ++i) {
      m.insert(rows.at(i), cols.at(i)) = count.at(i);
    }
  }

  m.makeCompressed();
  return m;
}

template <typename N>
inline void Cooler::fetch_df(COODataFrame<N>& buff, std::string_view range1,
                             std::string_view range2, std::string_view normalization) {
  if (!_clr) {
    throw std::runtime_error("Cooler::fetch_df() was called on an un-initialized object");
  }

  const auto divisive_weights =
      infer_weight_type(uri(), normalization) == balancing::Weights::Type::DIVISIVE;

  auto selector =
      normalization == "NONE"
          ? _clr.attr("matrix")("count", false, false, true, false, true, divisive_weights)
          : _clr.attr("matrix")("count", normalization, false, true, false, true, divisive_weights);

  if (range2.empty()) {
    range2 = range1;
  }

  buff = selector.attr("fetch")(range1, range2);
}

template <typename N>
inline void Cooler::fetch_df(BG2DataFrame<N>& buff, std::string_view range1,
                             std::string_view range2, std::string_view normalization) {
  if (!_clr) {
    throw std::runtime_error("Cooler::fetch_df() was called on an un-initialized object");
  }

  const auto divisive_weights =
      infer_weight_type(uri(), normalization) == balancing::Weights::Type::DIVISIVE;

  auto selector =
      normalization == "NONE"
          ? _clr.attr("matrix")("count", false, false, true, true, true, divisive_weights)
          : _clr.attr("matrix")("count", normalization, false, true, true, true, divisive_weights);

  if (range2.empty()) {
    range2 = range1;
  }

  buff = selector.attr("fetch")(range1, range2);
}

template <typename N>
inline Eigen2DDense<N> Cooler::fetch_dense(std::string_view range1, std::string_view range2,
                                           std::string_view normalization) {
  if (!_clr) {
    throw std::runtime_error("Cooler::fetch_dense() was called on an un-initialized object");
  }

  if (normalization != "NONE" && std::is_integral_v<N>) {
    throw std::runtime_error(
        "fetching balanced interactions requires Eigen2DDense<N> to be of floating-point type");
  }

  const auto divisive_weights =
      infer_weight_type(uri(), normalization) == balancing::Weights::Type::DIVISIVE;

  auto selector =
      normalization == "NONE"
          ? _clr.attr("matrix")("count", false, false, false, false, true, divisive_weights)
          : _clr.attr("matrix")("count", normalization, false, false, false, true,
                                divisive_weights);

  if (range2.empty()) {
    range2 = range1;
  }

  using NumpyArray = pybind11::array_t<N, pybind11::array::c_style | pybind11::array::forcecast>;
  const auto m_py = pybind11::cast<NumpyArray>(selector.attr("fetch")(range1, range2));

  const auto src = m_py.template unchecked<2>();
  Eigen2DDense<N> dest(m_py.shape(0), m_py.shape(1));
  for (std::int64_t i = 0; i < dest.rows(); ++i) {
    for (std::int64_t j = 0; j < dest.cols(); ++j) {
      dest(i, j) = src(i, j);
    }
  }

  return dest;
}

template <typename N>
inline EigenSparse<N> Cooler::fetch_sparse(std::string_view range1, std::string_view range2,
                                           std::string_view normalization) {
  if (!_clr) {
    throw std::runtime_error("Cooler::fetch_sparse() was called on an un-initialized object");
  }

  if (normalization != "NONE" && std::is_integral_v<N>) {
    throw std::runtime_error(
        "fetching balanced interactions requires EigenSparse<N> to be of floating-point type");
  }

  const auto divisive_weights =
      infer_weight_type(uri(), normalization) == balancing::Weights::Type::DIVISIVE;

  auto selector =
      normalization == "NONE"
          ? _clr.attr("matrix")("count", false, true, false, false, true, divisive_weights)
          : _clr.attr("matrix")("count", normalization, true, false, false, true, divisive_weights);

  if (range2.empty()) {
    range2 = range1;
  }

  return scipy_coo_to_eigen<N>(selector.attr("fetch")(range1, range2));
}

inline balancing::Weights::Type Cooler::infer_weight_type(std::string_view uri,
                                                          std::string_view normalization) {
  if (normalization == "NONE") {
    return balancing::Weights::Type::MULTIPLICATIVE;
  }
  const auto weights = hictk::cooler::File{uri}.normalization(normalization);
  if (!weights) {
    return balancing::Weights::Type::UNKNOWN;
  }

  return weights->type();
}

}  // namespace hictk::fuzzer::cooler
