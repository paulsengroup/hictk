// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/transformers/coarsen.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

using namespace hictk::transformers;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (cooler): coarsen", "[transformers][short]") {
  SECTION("simple") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/200000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch("1");
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = clr2.fetch("1");

    const auto v1 = sel1.read_all();
    const auto v2 = sel2.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  SECTION("coarsen recursive") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/400000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch("1");
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = CoarsenPixels(sel1.begin(), sel1.end(), sel1.dest_bins_ptr(), 2);
    auto sel3 = clr2.fetch("1");

    const auto v1 = sel2.read_all();
    const auto v2 = sel3.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  SECTION("coarsen gw") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/200000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch();
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = clr2.fetch();

    const auto v1 = sel1.read_all();
    const auto v2 = sel2.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  SECTION("coarsen empty range") {
    const auto path = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const cooler::File clr1(path.string());

    auto sel = clr1.fetch();
    auto sel1 = CoarsenPixels(sel.end<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);

    CHECK(sel1.begin() == sel1.end());
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (hic): coarsen", "[transformers][short]") {
  auto path = datadir / "hic/4DNFIZ1ZVXC8.hic8";

  const hic::File hf1(path.string(), 500'000);
  const hic::File hf2(path.string(), 2'500'000);

  auto sel = hf1.fetch("chr2R");
  auto sel1 = CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), hf1.bins_ptr(), 5);
  auto sel2 = hf2.fetch("chr2R");

  const auto v1 = sel1.read_all();
  const auto v2 = sel2.read_all<std::int32_t>();
  REQUIRE(v1.size() == v2.size());

  for (std::size_t i = 0; i < v1.size(); ++i) {
    CHECK(v1[i] == v2[i].to_thin());
  }
}

}  // namespace hictk::test::transformers
