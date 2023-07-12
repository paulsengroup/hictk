// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/transformers.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>

#include "hictk/cooler.hpp"
#include "hictk/hic.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

using namespace hictk::transformers;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (cooler)", "[transformers][short]") {
  auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
  auto clr = cooler::File::open_read_only(path.string());

  SECTION("join genomic coords") {
    auto sel = clr.fetch("chr1", 5'000'000, 10'000'000);
    const auto jsel = transformers::JoinGenomicCoords(sel.begin<std::int32_t>(),
                                                      sel.end<std::int32_t>(), clr.bins_ptr());
    constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
    const auto pixels = jsel.read_all();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(pixels[i].coords.bin1.start() == expected[i]);
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (hic)", "[transformers][short]") {
  auto path = datadir / "hic/4DNFIZ1ZVXC8.hic8";
  auto hf = hic::HiCFile(path.string(), 2'500'000);

  SECTION("join genomic coords") {
    auto sel = hf.fetch("chr2L", 5'000'000, 10'000'000);
    const auto jsel = transformers::JoinGenomicCoords(sel.begin<std::int32_t>(),
                                                      sel.end<std::int32_t>(), hf.bins_ptr());
    constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
    const auto pixels = jsel.read_all();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(pixels[i].coords.bin1.start() == expected[i]);
    }
  }
}

}  // namespace hictk::test::transformers
