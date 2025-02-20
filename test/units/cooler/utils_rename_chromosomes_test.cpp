// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <map>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: utils rename chromosomes", "[cooler][short]") {
  const Reference ref2{{0, "1", 10}, {1, "2", 10}};
  const Reference ref3{{0, "chr1", 10}, {1, "2", 10}};

  const auto path = testdir() / "rename_chromosomes.cool";

  const Reference ref{{0, "chr1", 10}, {1, "chr2", 10}};
  std::ignore = cooler::File::create(path.string(), ref, 1, true);
  cooler::utils::rename_chromosomes(path.string(),
                                    std::map<std::string, std::string>{{"chr1", "1"}});

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

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::cooler_file
