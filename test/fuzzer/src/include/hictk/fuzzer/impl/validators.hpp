// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

#include "hictk/fuzzer/cooler.hpp"
#include "hictk/pixel.hpp"

namespace hictk::fuzzer {

namespace internal {

template <typename N>
constexpr bool is_close(N n1, N n2, double rtol, double atol) noexcept {
  assert(rtol >= 0 && rtol <= 1);
  if constexpr (std::is_integral_v<N>) {
    return n1 == n2;
  } else {
    if (n1 == n2) {
      return true;
    }
    if (std::isnan(n1)) {
      return std::isnan(n2);
    }

    // https://peps.python.org/pep-0485/
    const auto diff = std::abs(n1 - n2);
    return (diff <= abs(rtol * n2)) || (diff <= atol);
  }
}

}  // namespace internal

template <typename N>
inline bool compare_pixels([[maybe_unused]] std::uint16_t task_id, std::string_view range1,
                           std::string_view range2, const std::vector<ThinPixel<N>>& expected,
                           const std::vector<ThinPixel<N>>& found) {
  if (expected.size() != found.size()) {
    SPDLOG_WARN(FMT_STRING("[{}]: {}, {}: FAIL! Expected {} nnz, found {}!"), task_id, range1,
                range2, expected.size(), found.size());
    return false;
  }

  std::size_t num_mismatches{};
  for (std::size_t i = 0; i < expected.size(); ++i) {
    const auto& p1 = expected[i];
    const auto& p2 = found[i];
    num_mismatches += p1.bin1_id != p2.bin1_id || p1.bin2_id != p2.bin2_id ||
                      !internal::is_close(p1.count, p2.count);
  }

  if (num_mismatches != 0) {
    SPDLOG_WARN(FMT_STRING("[{}]: {}, {}: FAIL! Found {} differences!"), task_id, range1, range2,
                num_mismatches);
    return false;
  }
  return true;
}

template <typename N>
inline bool compare_pixels([[maybe_unused]] std::uint16_t task_id, std::string_view range1,
                           std::string_view range2, const std::vector<Pixel<N>>& expected,
                           const std::vector<Pixel<N>>& found) {
  if (expected.size() != found.size()) {
    SPDLOG_WARN(FMT_STRING("[{}]: {}, {}: FAIL! Expected {} nnz, found {}!"), task_id, range1,
                range2, expected.size(), found.size());
    return false;
  }

  std::size_t num_mismatches{};
  for (std::size_t i = 0; i < expected.size(); ++i) {
    const auto& p1 = expected[i];
    const auto& p2 = found[i];
    num_mismatches += p1.coords.bin1 != p2.coords.bin1 || p1.coords.bin2 != p2.coords.bin2 ||
                      !internal::is_close(p1.count, p2.count);
  }

  if (num_mismatches != 0) {
    SPDLOG_WARN(FMT_STRING("[{}]: {}, {}: FAIL! Found {} differences!"), task_id, range1, range2,
                num_mismatches);
    return false;
  }
  return true;
}

template <typename N>
inline bool compare_pixels([[maybe_unused]] std::uint16_t task_id, std::string_view range1,
                           std::string_view range2, const Eigen2DDense<N>& expected,
                           const Eigen2DDense<N>& found) {
  if (expected.rows() != found.rows() || expected.cols() != found.cols()) {
    SPDLOG_WARN(
        FMT_STRING("[{}]: {}, {}: FAIL! Expected matrix of shape [{}, {}], found [{}, {}]!"),
        task_id, range1, range2, expected.rows(), expected.cols(), found.rows(), found.cols());
    return false;
  }

  std::size_t num_mismatches{};
  for (std::int64_t i = 0; i < expected.size(); ++i) {
    num_mismatches += !internal::is_close(*(expected.data() + i), *(found.data() + i));
  }

  if (num_mismatches != 0) {
    SPDLOG_WARN(FMT_STRING("[{}]: {}, {}: FAIL! Found {} differences!"), task_id, range1, range2,
                num_mismatches);
    return false;
  }
  return true;
}

template <typename N>
inline bool compare_pixels([[maybe_unused]] std::uint16_t task_id, std::string_view range1,
                           std::string_view range2, const EigenSparse<N>& expected,
                           const EigenSparse<N>& found) {
  if (expected.rows() != found.rows() || expected.cols() != found.cols()) {
    SPDLOG_WARN(
        FMT_STRING("[{}]: {}, {}: FAIL! Expected matrix of shape [{}, {}], found [{}, {}]!"),
        task_id, range1, range2, expected.rows(), expected.cols(), found.rows(), found.cols());
    return false;
  }

  if (expected.nonZeros() != found.nonZeros()) {
    SPDLOG_WARN(FMT_STRING("[{}]: {}, {}: FAIL! Expected {} nnz, found {}!"), task_id, range1,
                range2, expected.nonZeros(), found.nonZeros());
    return false;
  }

  // FIXME this doesn't work because cooler mirrors interactions even when returning them as sparse
  // matrices
  return compare_pixels(task_id, range1, range2, Eigen2DDense<N>{expected.toDense()},
                        Eigen2DDense<N>{found.toDense()});
}

}  // namespace hictk::fuzzer
