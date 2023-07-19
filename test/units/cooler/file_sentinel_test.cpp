// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "hictk/cooler/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: sentinel attribute", "[cooler][short]") {
  const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};

  const auto path = testdir() / "test_sentinel_attr.cool";
  constexpr std::uint32_t bin_size = 1000;
  auto f = File::create(path.string(), chroms, bin_size, true);

  SECTION("Read-only") {
    const auto path1 = datadir / "cooler_test_file.cool";
    const auto f1 = File::open(path1.string());
    CHECK(Attribute::read<std::uint8_t>(f1.group("/")(), internal::SENTINEL_ATTR_NAME) !=
          internal::SENTINEL_ATTR_VALUE);
  }

  SECTION("Create") {
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) ==
          internal::SENTINEL_ATTR_VALUE);
    f.close();
    f = File::open(path.string());
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) !=
          internal::SENTINEL_ATTR_VALUE);
  }

  SECTION("Create (file was not closed properly)") {
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) ==
          internal::SENTINEL_ATTR_VALUE);

    CHECK_THROWS(f = File::open(path.string()));
    CHECK_THROWS(f = File::create(path.string(), chroms, bin_size, true));
  }
}

}  // namespace hictk::cooler::test::cooler_file
