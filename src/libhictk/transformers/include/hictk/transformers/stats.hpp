// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <type_traits>

namespace hictk::transformers {

template <typename PixelIt>
[[nodiscard]] double avg(PixelIt first, PixelIt last);

template <typename PixelIt, typename N = std::decay_t<decltype(std::declval<PixelIt>()->count)>>
[[nodiscard]] N max(PixelIt first, PixelIt last);

template <typename PixelIt>
[[nodiscard]] std::size_t nnz(PixelIt first, PixelIt last);

template <typename PixelIt, typename N = std::decay_t<decltype(std::declval<PixelIt>()->count)>>
[[nodiscard]] N sum(PixelIt first, PixelIt last);

}  // namespace hictk::transformers

#include "./impl/stats_impl.hpp"  // NOLINT
