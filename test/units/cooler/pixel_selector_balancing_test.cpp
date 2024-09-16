// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>

#include "hictk/balancing/weights.hpp"
#include "hictk/cooler/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::pixel_selector {
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: pixel selector w/ balancing", "[pixel_selector][short]") {
  auto path = datadir / "ENCFF993FGR.2500000.cool";
  File clr(path.string());

  SECTION("read weights") {
    SECTION("valid") {
      CHECK(clr.normalization("weight").type() == hictk::balancing::Weights::Type::MULTIPLICATIVE);
      for (const auto* name : {"GW_SCALE", "INTER_SCALE", "SCALE", "VC", "VC_SQRT"}) {
        CHECK(clr.normalization(name).type() == hictk::balancing::Weights::Type::DIVISIVE);
      }
    }

    SECTION("invalid") {
      CHECK_THROWS(clr.normalization(""));
      CHECK_THROWS(clr.normalization("AAA"));
    }

    SECTION("purging") {
      CHECK(clr.purge_weights() == false);
      CHECK(clr.purge_weights("weight") == false);

      const auto w = clr.normalization_ptr("weight");
      CHECK(w.use_count() == 2);
      CHECK(clr.purge_weights("weight") == true);
      CHECK(w.use_count() == 1);

      std::ignore = clr.normalization("weight");
      CHECK(clr.purge_weights() == true);
    }
  }

  SECTION("1D query") {
    const auto selector = clr.fetch("chr1", 5'000'000, 10'000'000, clr.normalization_ptr("weight"));
    constexpr std::array<double, 3> expected{3.345797, 0.328794, 4.456354};
    const auto pixels = selector.read_all<double>();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK_THAT(pixels[i].count, Catch::Matchers::WithinAbs(expected[i], 1.0e-6));
    }
  }

  SECTION("2D query") {
    const auto selector = clr.fetch("chr1", 5'000'000, 10'000'000, "chr2", 5'000'000, 10'000'000,
                                    clr.normalization_ptr("weight"));
    constexpr std::array<double, 4> expected{0.001782, 0.002756, 0.002047, 0.004749};
    const auto pixels = selector.read_all<double>();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK_THAT(pixels[i].count, Catch::Matchers::WithinAbs(expected[i], 1.0e-6));
    }
  }

  SECTION("invalid iterator type") {
    const auto selector = clr.fetch("chr1", 5'000'000, 10'000'000, "chr2", 5'000'000, 10'000'000,
                                    clr.normalization_ptr("weight"));
    CHECK_THROWS(selector.read_all<std::int32_t>());
  }
}

}  // namespace hictk::cooler::test::pixel_selector
