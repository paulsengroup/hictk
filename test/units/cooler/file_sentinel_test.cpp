// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>

#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler.hpp"
#include "hictk/reference.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::cooler_file {

static const auto& datadir = hictk::test::datadir;
static const auto& testdir = hictk::test::testdir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: sentinel attribute", "[cooler][short]") {
  const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};

  const auto path = testdir() / "test_sentinel_attr.cool";
  constexpr std::uint32_t bin_size = 1000;
  auto f = File::create(path.string(), chroms, bin_size, true);

  SECTION("Read-only") {
    const auto path1 = datadir / "cooler" / "cooler_test_file.cool";
    const File f1(path1.string());
    CHECK(Attribute::read<std::uint8_t>(f1.group("/")(), internal::SENTINEL_ATTR_NAME) !=
          internal::SENTINEL_ATTR_VALUE);
  }

  SECTION("Create") {
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) ==
          internal::SENTINEL_ATTR_VALUE);
    f.close();
    f = File(path.string());
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) !=
          internal::SENTINEL_ATTR_VALUE);
  }

  SECTION("Create (file was not closed properly)") {
    CHECK(Attribute::read<std::uint8_t>(f.group("/")(), internal::SENTINEL_ATTR_NAME) ==
          internal::SENTINEL_ATTR_VALUE);

    CHECK_THROWS(f = File(path.string()));
    CHECK_THROWS(f = File::create(path.string(), chroms, bin_size, true));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::cooler_file
