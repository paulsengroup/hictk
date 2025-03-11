// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "hictk/cooler/utils.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::utils {
static const auto& testdir = hictk::test::testdir;
static const auto& datadir = hictk::test::datadir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: utils equal", "[equal][utils][short]") {
  const auto path1 = datadir / "cooler" / "cooler_test_file.cool";
  const auto path2 = datadir / "cooler" / "multires_cooler_test_file.mcool::/resolutions/6400000";

  const auto path3 = testdir() / "cooler_equal_test.cool";

  std::filesystem::remove(path3);
  std::filesystem::copy(path1, path3);

  SECTION("equal") {
    CHECK(cooler::utils::equal(path1.string(), path1.string()));
    CHECK(cooler::utils::equal(path1.string(), path3.string()));
    CHECK(cooler::utils::equal(path1.string(), path3.string(), false));
  }

  SECTION("not equal") { CHECK_FALSE(cooler::utils::equal(path1.string(), path2.string())); }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::utils
