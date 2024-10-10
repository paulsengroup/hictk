// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>

#include "hictk/cooler/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::pixel_selector {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: pixel selector 2D queries", "[pixel_selector][short]") {
  using T = std::uint32_t;
  const auto path = datadir / "cooler_test_file.cool";
  const File f(path.string());

  SECTION("cis") {
    SECTION("overloads return identical results") {
      const auto sel1 = f.fetch("1", "1");
      const auto sel2 = f.fetch("1", 0, 197195432, "1", 0, 197195432);
      CHECK(sel1 == sel2);

      const auto& pixels1 = sel1.read_all<T>();
      const auto& pixels2 = sel2.read_all<T>();

      REQUIRE(pixels1.size() == pixels2.size());
    }

    SECTION("valid") {
      auto selector = f.fetch("1:5000000-5500000", "1:5000000-6500000");
      const auto pixels = selector.read_all<T>();
      REQUIRE(pixels.size() == 8);

      CHECK(pixels[0].count == 20);
      CHECK(pixels[1].count == 1);
      CHECK(pixels[2].count == 18);
      CHECK(pixels[3].count == 8);
      CHECK(pixels[4].count == 1);
      CHECK(pixels[5].count == 9);
      CHECK(pixels[6].count == 6);
      CHECK(pixels[7].count == 2);
    }

    SECTION("invalid") {
      SECTION("query overlaps lower triangle") {
        CHECK_THROWS_WITH(f.fetch("1:6000000-6500000", "1:5000000-5500000"),
                          Catch::Matchers::ContainsSubstring("overlaps with the lower-triangle"));
      }
    }

    SECTION("empty") {
      auto selector = f.fetch("1:0-100000");
      CHECK(selector.begin<T>() == selector.end<T>());
    }
  }

  SECTION("trans") {
    SECTION("overloads return identical results") {
      const auto sel1 = f.fetch("1", "4");
      const auto sel2 = f.fetch("1", 0, 197195432, "4", 0, 155630120, nullptr);

      CHECK(sel1 == sel2);

      const auto& pixels1 = sel1.read_all<T>();
      const auto& pixels2 = sel2.read_all<T>();

      REQUIRE(pixels1.size() == pixels2.size());
    }
    SECTION("valid") {
      auto selector = f.fetch("1:48000000-50000000", "4:30000000-35000000");
      const auto pixels = selector.read_all<T>();
      REQUIRE(pixels.size() == 6);

      CHECK(pixels[0].count == 1);
      CHECK(pixels[1].count == 3);
      CHECK(pixels[2].count == 1);
      CHECK(pixels[3].count == 3);
      CHECK(pixels[4].count == 7);
      CHECK(pixels[5].count == 1);
    }

    SECTION("invalid") {
      SECTION("query overlaps lower triangle") {
        CHECK_THROWS_WITH(f.fetch("2", "1"),
                          Catch::Matchers::ContainsSubstring("overlaps with the lower-triangle"));
      }
    }

    SECTION("empty") {
      auto selector = f.fetch("1:0-50000", "2:0-50000");
      CHECK(selector.begin<T>() == selector.end<T>());
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::pixel_selector
