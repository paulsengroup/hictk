// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/reference.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: read/write bin table", "[cooler][short]") {
  const auto path = (testdir() / "test_write_bin_table.cool").string();

  const Reference chroms{Chromosome{0, "chr1", 50001}, Chromosome{1, "chr2", 25017},
                         Chromosome{2, "chr3", 10000}};

  constexpr std::uint32_t bin_size = 5000;
  const BinTable table(chroms, bin_size);

  {
    auto f = File::create(path, chroms, bin_size, true);
  }

  File f(path);

  auto start_it = f.dataset("bins/start").begin<std::uint32_t>(32'000);
  auto end_it = f.dataset("bins/end").begin<std::uint32_t>(32'000);

  REQUIRE(start_it != f.dataset("bins/start").end<std::uint32_t>(0));
  REQUIRE(end_it != f.dataset("bins/end").end<std::uint32_t>(0));

  for (const auto bin : table) {
    CHECK(*start_it++ == bin.start());
    CHECK(*end_it++ == bin.end());
  }

  CHECK(start_it == f.dataset("bins/start").end<std::uint32_t>(0));
  CHECK(end_it == f.dataset("bins/end").end<std::uint32_t>(0));
}

TEST_CASE("Cooler: validate bin table", "[cooler][short]") {
  auto path = datadir / "ENCFF993FGR.2500000.cool";

  const File f(path.string());

  CHECK_NOTHROW(f.validate_bins(true));
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::cooler_file
