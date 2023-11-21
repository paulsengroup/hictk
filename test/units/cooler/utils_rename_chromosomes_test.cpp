// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: rename chromosomes", "[cooler][short]") {
  const Reference ref2{{0, "1", 10}, {1, "2", 10}};
  const Reference ref3{{0, "chr1", 10}, {1, "2", 10}};

  const auto path = testdir() / "rename_chromosomes.cool";

  const Reference ref{{0, "chr1", 10}, {1, "chr2", 10}};
  std::ignore = cooler::File::create(path.string(), ref, 1, true);
  cooler::utils::rename_chromosomes(path.string(), {{"chr1", "1"}});

  {
    const auto chroms = cooler::File(path.string()).chromosomes();
    CHECK(chroms.size() == 2);
    CHECK(!chroms.contains("chr1"));
    CHECK(chroms.contains("chr2"));
    CHECK(chroms.contains("1"));
  }

  const std::vector<std::pair<std::string, std::string>> mappings{{"1", "abc12345"}};
  cooler::utils::rename_chromosomes(path.string(), mappings.begin(), mappings.end());

  {
    const auto chroms = cooler::File(path.string()).chromosomes();
    CHECK(chroms.size() == 2);
    CHECK(!chroms.contains("1"));
    CHECK(chroms.contains("chr2"));
    CHECK(chroms.contains("abc12345"));
  }
}

}  // namespace hictk::cooler::test::cooler_file
