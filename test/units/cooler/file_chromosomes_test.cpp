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
TEST_CASE("Cooler: read/write chromosomes", "[cooler][short]") {
  const auto path = (testdir() / "test_write_chroms.cool").string();

  constexpr std::uint32_t bin_size = 5000;
  const Reference chroms{Chromosome{0, "chr1", 50001}, Chromosome{1, "chr2", 25017},
                         Chromosome{2, "chr3", 10000}};

  {
    auto f = File::create(path, chroms, bin_size, true);
    CHECK(chroms == f.chromosomes());
  }

  const auto f = File::open(path, DEFAULT_HDF5_CACHE_SIZE, false);
  CHECK(chroms == f.chromosomes());
}

}  // namespace hictk::cooler::test::cooler_file
