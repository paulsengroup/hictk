// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <numeric>

namespace hictk::transformers {

template <typename PixelIt, typename N>
inline N sum(PixelIt first, PixelIt last) {
  return std::accumulate(first, last, N(0),
                         [&](const N accumulator, const auto& p) { return accumulator + p.count; });
}

template <typename PixelIt, typename N>
inline N max(PixelIt first, PixelIt last) {
  // I prefer using for_each instead of max_element to avoid keeping around a copy
  // of the iterator for the current max_element
  N max_ = 0;
  std::for_each(first, last, [&](const auto& p) { max_ = std::max(max_, p.count); });
  return max_;
}

template <typename PixelIt>
inline std::size_t nnz(PixelIt first, PixelIt last) {
  return static_cast<std::size_t>(std::distance(first, last));
}

template <typename PixelIt>
inline double avg(PixelIt first, PixelIt last) {
  std::size_t nnz = 0;
  const auto sum =
      std::accumulate(first, last, double{0}, [&](const double accumulator, const auto& p) {
        ++nnz;
        return accumulator + p.count;
      });
  if (nnz == 0) {
    return 0.0;
  }
  return sum / static_cast<double>(nnz);
}

}  // namespace hictk::transformers
