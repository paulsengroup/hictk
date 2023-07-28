// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>

#include "hictk/cooler/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::pixel_selector {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: pixel selector 2D queries", "[pixel_selector][short]") {
  using T = std::uint32_t;
  const auto path = datadir / "cooler_test_file.cool";
  const File f(path.string());

  SECTION("cis") {
    SECTION("overloads return identical results") {
      CHECK(f.fetch("1:5000000-5500000", "1:5000000-6500000") ==
            f.fetch("1", 5000000, 5500000, "1", 5000000, 6500000));
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

    SECTION("empty") {
      auto selector = f.fetch("1:0-100000");
      CHECK(selector.begin<T>() == selector.end<T>());
    }
  }

  SECTION("trans") {
    SECTION("overloads return identical results") {
      CHECK(f.fetch("1:48000000-50000000", "4:30000000-35000000") ==
            f.fetch("1", 48000000, 50000000, "4", 30000000, 35000000));
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

    SECTION("empty") {
      auto selector = f.fetch("1:0-50000", "2:0-50000");
      CHECK(selector.begin<T>() == selector.end<T>());
    }
  }
}

}  // namespace hictk::cooler::test::pixel_selector
