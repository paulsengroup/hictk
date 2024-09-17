// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <cstdint>
#include <vector>

#include "hictk/balancing/weights.hpp"

namespace hictk::test::balancing {
inline void compare_weights(const hictk::balancing::Weights& weights_,
                            const hictk::balancing::Weights& expected_, double atol = 1.0e-5,
                            double rtol = 1.0e-5) {
  REQUIRE(weights_.size() == expected_.size());

  const auto weights = weights_(hictk::balancing::Weights::Type::DIVISIVE);
  const auto expected = expected_(hictk::balancing::Weights::Type::DIVISIVE);

  for (std::size_t i = 0; i < weights.size(); ++i) {
    if (std::isnan(expected[i])) {
      CHECK(std::isnan(weights[i]));
    } else {
      // Basically we don't care about the relative error when the weights are very small, as this
      // will not lead to significant differences when balancing interactions
      const auto delta = std::abs(weights[i] - expected[i]);
      if (delta > atol) {
        CHECK_THAT(weights[i], Catch::Matchers::WithinRel(expected[i], rtol));
      } else {
        CHECK_THAT(weights[i], Catch::Matchers::WithinAbs(expected[i], atol));
      }
    }
  }
}

template <typename T>
inline void compare_vectors(const std::vector<T>& v1, const std::vector<T>& v2) {
  REQUIRE(v1.size() == v2.size());

  for (std::size_t i = 0; i < v1.size(); ++i) {
    CHECK(v1[i] == v2[i]);
  }
}

}  // namespace hictk::test::balancing
