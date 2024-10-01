// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/transformers/join_genomic_coords.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

using namespace hictk::transformers;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (cooler): join genomic coords", "[transformers][short]") {
  const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
  const cooler::File clr(path.string());

  auto sel = clr.fetch("chr1", 5'000'000, 10'000'000);
  SECTION("range with data") {
    const auto jsel =
        JoinGenomicCoords(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr.bins_ptr());
    constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
    const auto pixels = jsel.read_all();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(pixels[i].coords.bin1.start() == expected[i]);
    }
  }
  SECTION("empty range") {
    const auto jsel =
        JoinGenomicCoords(sel.end<std::int32_t>(), sel.end<std::int32_t>(), clr.bins_ptr());
    CHECK(jsel.begin() == jsel.end());
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (hic): join genomic coords", "[transformers][short]") {
  auto path = datadir / "hic/4DNFIZ1ZVXC8.hic8";

  const hic::File hf(path.string(), 2'500'000);
  auto sel = hf.fetch("chr2L", 5'000'000, 10'000'000);
  const auto jsel =
      JoinGenomicCoords(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), hf.bins_ptr());
  constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
  const auto pixels = jsel.read_all();
  REQUIRE(pixels.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    CHECK(pixels[i].coords.bin1.start() == expected[i]);
  }
}

}  // namespace hictk::test::transformers
