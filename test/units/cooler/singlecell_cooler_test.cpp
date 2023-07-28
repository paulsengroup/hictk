// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/singlecell_cooler.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "tmpdir.hpp"

namespace hictk::cooler::test::singlecell_cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("SingleCellCooler: open read-only", "[cooler][short]") {
  const auto path = datadir / "single_cell_cooler_test_file.scool";

  const SingleCellFile sclr(path.string());

  CHECK(sclr.attributes().format == SCOOL_MAGIC);
  CHECK(sclr.attributes().format_version == 1);
  CHECK(sclr.attributes().ncells == 5);

  const auto& cells = sclr.cells();
  REQUIRE(cells.size() == 5);

  const auto first_cell_name = *cells.begin();
  CHECK(utils::is_cooler(sclr.open(first_cell_name).uri()));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("SingleCellCooler: create cells", "[cooler][short]") {
  const auto base_path = datadir / "cooler_test_file.cool";
  const File base_clr(base_path.string());

  const auto path = testdir() / "test_create_cells.scool";

  auto sclr =
      SingleCellFile::create(path.string(), base_clr.chromosomes(), base_clr.bin_size(), true);

  SECTION("valid cells") {
    for (std::size_t i = 0; i < 10; ++i) {
      const auto uri = sclr.create_cell<std::int32_t>(std::to_string(i)).uri();
      CHECK(utils::is_cooler(uri));
    }

    CHECK(sclr.cells().size() == 10);
  }

  SECTION("duplicate cell") {
    sclr.create_cell<std::int32_t>("A");
    CHECK_THROWS(sclr.create_cell<std::int32_t>("A"));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("SingleCellCooler: aggregate cells", "[cooler][short]") {
  const auto base_path = datadir / "cooler_test_file.cool";
  const File base_clr(base_path.string());

  const auto path1 = testdir() / "test_aggregate_cells.scool";
  const auto path2 = testdir() / "test_aggregate_cells.cool";

  auto sclr =
      SingleCellFile::create(path1.string(), base_clr.chromosomes(), base_clr.bin_size(), true);

  {
    auto clr1 = sclr.create_cell<std::int32_t>("A");
    auto clr2 = sclr.create_cell<std::int32_t>("B");

    clr1.append_pixels(base_clr.begin<std::int32_t>(), base_clr.end<std::int32_t>());
    clr2.append_pixels(base_clr.begin<std::int32_t>(), base_clr.end<std::int32_t>());
  }

  sclr.aggregate<std::int32_t>(path2.string());

  const File clr(path2.string());
  CHECK(std::get<std::int64_t>(*clr.attributes().sum) ==
        2 * std::get<std::int64_t>(*base_clr.attributes().sum));
}

}  // namespace hictk::cooler::test::singlecell_cooler_file
