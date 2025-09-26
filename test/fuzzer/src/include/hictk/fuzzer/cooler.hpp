// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <cstdint>
#include <memory>
#include <string_view>
#include <type_traits>
#include <vector>

#include "hictk/balancing/weights.hpp"
#include "hictk/fuzzer/common.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk::fuzzer::cooler {

[[nodiscard]] std::string_view version();

template <typename N>
struct COODataFrame {
  static_assert(std::is_arithmetic_v<N>);
  NumpyArray<std::int64_t> bin1_id{};
  NumpyArray<std::int64_t> bin2_id{};

  NumpyArray<N> count{};

  // NOLINTNEXTLINE(*-unnecessary-value-param)
  COODataFrame<N>& operator=(pybind11::object df);

  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] std::vector<ThinPixel<N>> to_vector() const;
  void to_vector(std::vector<ThinPixel<N>>& buffer) const;
};

template <typename N>
struct BG2DataFrame {
  static_assert(std::is_arithmetic_v<N>);
  pybind11::list chrom1{};
  NumpyArray<std::int32_t> start1{};
  NumpyArray<std::int32_t> end1{};

  pybind11::list chrom2{};
  NumpyArray<std::int32_t> start2{};
  NumpyArray<std::int32_t> end2{};

  NumpyArray<N> count{};

  // NOLINTNEXTLINE(*-unnecessary-value-param)
  BG2DataFrame<N>& operator=(pybind11::object df);

  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] std::vector<Pixel<N>> to_vector(const Reference& chroms) const;
  void to_vector(const Reference& chroms, std::vector<Pixel<N>>& buffer) const;
};

template <typename N>  // NOLINTNEXTLINE(*-unnecessary-value-param)
[[nodiscard]] EigenSparse<N> scipy_coo_to_eigen(pybind11::object obj);

class Cooler {
  pybind11::object _clr{};

 public:
  Cooler() = default;
  explicit Cooler(std::string_view uri);

  [[nodiscard]] std::string uri() const;
  [[nodiscard]] std::uint32_t resolution() const;

  template <typename N>
  void fetch_df(COODataFrame<N>& buff, std::string_view range1, std::string_view range2 = "",
                std::string_view normalization = "NONE",
                std::optional<std::uint64_t> diagonal_band_width = {});
  template <typename N>
  void fetch_df(BG2DataFrame<N>& buff, std::string_view range1, std::string_view range2 = "",
                std::string_view normalization = "NONE",
                std::optional<std::uint64_t> diagonal_band_width = {});
  template <typename N>
  [[nodiscard]] Eigen2DDense<N> fetch_dense(std::string_view range1, std::string_view range2 = "",
                                            std::string_view normalization = "NONE");
  template <typename N>
  [[nodiscard]] EigenSparse<N> fetch_sparse(std::string_view range1, std::string_view range2,
                                            std::string_view normalization);

  [[nodiscard]] static balancing::Weights::Type infer_weight_type(std::string_view uri,
                                                                  std::string_view normalization);
};

}  // namespace hictk::fuzzer::cooler

#include "./impl/cooler.hpp"
