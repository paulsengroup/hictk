// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>
#include <vector>

namespace hictk::fuzzer {

namespace internal {

template <typename N>
[[nodiscard]] constexpr bool is_close(N n1, N n2, double rtol = 1.0e-6, double atol = 0) noexcept;

}  // namespace internal

template <typename N>
[[nodiscard]] bool compare_pixels(std::string_view range1, std::string_view range2,
                                  const std::vector<ThinPixel<N>>& expected,
                                  const std::vector<ThinPixel<N>>& found);

template <typename N>
[[nodiscard]] bool compare_pixels(std::string_view range1, std::string_view range2,
                                  const std::vector<Pixel<N>>& expected,
                                  const std::vector<Pixel<N>>& found);
}  // namespace hictk::fuzzer

#include "./impl/validators.hpp"
