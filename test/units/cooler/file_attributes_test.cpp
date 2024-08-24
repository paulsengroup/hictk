// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <variant>

#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: read attributes", "[cooler][short]") {
  auto path = datadir / "cooler_test_file.cool";
  const File f(path.string());

  SECTION("bin size") { CHECK(f.resolution() == 100'000); }

  SECTION("common attributes") {
    const auto& attrs = f.attributes();
    CHECK(attrs.bin_size == 100'000);
    CHECK(attrs.bin_type == BinTable::Type::fixed);
    CHECK(attrs.creation_date == "2020-07-08T13:41:20.376258");
    CHECK(attrs.format == COOL_MAGIC);
    CHECK(attrs.format_url == "https://github.com/mirnylab/cooler");
    CHECK(attrs.format_version == 3);
    CHECK(attrs.generated_by == "cooler-0.8.8-dev");
    CHECK(attrs.assembly == "unknown");
    CHECK(attrs.metadata == "{}");
    CHECK(attrs.nbins == 26398);
    CHECK(attrs.nchroms == 20);
    CHECK(attrs.nnz == 107041);
    CHECK(attrs.storage_mode == "symmetric-upper");
    CHECK(attrs.sum.has_value());
    if (attrs.sum.has_value()) {
      std::visit([](auto& sum) { CHECK(sum == 395465); }, *attrs.sum);
    }
    CHECK(!attrs.cis.has_value());
  }
}

}  // namespace hictk::cooler::test::cooler_file
