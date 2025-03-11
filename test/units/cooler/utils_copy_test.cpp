// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::utils {
static const auto& testdir = hictk::test::testdir;
static const auto& datadir = hictk::test::datadir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: utils copy", "[copy][utils][short]") {
  SECTION("cooler -> cooler") {
    const auto src = datadir / "cooler" / "cooler_test_file.cool";
    const auto dest = testdir() / "cooler_copy_001.cool";

    cooler::utils::copy(src.string(), dest.string());
    CHECK(cooler::utils::equal(src.string(), dest.string()));
  }

  SECTION("cooler -> mcool") {
    const auto src = datadir / "cooler" / "cooler_test_file.cool";
    const auto dest = testdir() / "cooler_copy_002.mcool";
    const auto dest_uri = fmt::format(FMT_STRING("{}::/resolutions/1000"), dest.string());

    {
      auto mclr = MultiResFile::create(dest, File(src.string()).chromosomes(), true);
      mclr.init_resolution(1000);
    }

    cooler::utils::copy(src.string(), dest_uri);
    CHECK(cooler::utils::equal(src.string(), dest_uri));
  }

  SECTION("mcool -> cooler") {
    const auto src = datadir / "cooler" / "multires_cooler_test_file.mcool::/resolutions/100000";
    const auto dest = testdir() / "cooler_copy_003.cool";

    cooler::utils::copy(src.string(), dest.string());
    CHECK(cooler::utils::equal(src.string(), dest.string()));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::utils
