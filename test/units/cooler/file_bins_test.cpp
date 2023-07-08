// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "hictk/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: read/write bin table", "[cooler][short]") {
  const auto path = (testdir() / "test_write_bin_table.cool").string();

  const Reference chroms{Chromosome{0, "chr1", 50001}, Chromosome{1, "chr2", 25017},
                         Chromosome{2, "chr3", 10000}};

  constexpr std::uint32_t bin_size = 5000;
  const BinTable table(chroms, bin_size);

  { auto f = File::create_new_cooler(path, chroms, bin_size, true); }

  auto f = File::open_read_only(path);

  auto start_it = f.dataset("bins/start").begin<std::uint32_t>();
  auto end_it = f.dataset("bins/end").begin<std::uint32_t>();

  REQUIRE(start_it != f.dataset("bins/start").end<std::uint32_t>());
  REQUIRE(end_it != f.dataset("bins/end").end<std::uint32_t>());

  for (const auto bin : table) {
    CHECK(*start_it++ == bin.start());
    CHECK(*end_it++ == bin.end());
  }

  CHECK(start_it == f.dataset("bins/start").end<std::uint32_t>());
  CHECK(end_it == f.dataset("bins/end").end<std::uint32_t>());
}

}  // namespace hictk::cooler::test::cooler_file
