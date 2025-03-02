// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/pixel.hpp"
#include "hictk/transformers/diagonal_band.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)

using namespace hictk::transformers;

template <typename T, typename PixelSelector>
static std::vector<ThinPixel<T>> fetch_pixels(const PixelSelector& sel, std::uint64_t num_bins) {
  std::vector<ThinPixel<T>> buff;

  auto first = sel.template begin<T>();
  const auto last = sel.template end<T>();

  while (first != last) {
    if (first->bin2_id - first->bin1_id < num_bins) {
      buff.emplace_back(*first);
    }
    ++first;
  }

  return buff;
}

TEST_CASE("Transformers (cooler): diagonal band", "[transformers][short]") {
  const cooler::File clr((datadir / "cooler" / "cooler_test_file.cool").string());
  SECTION("simple") {
    constexpr std::uint64_t num_bins = 200;

    const auto sel = clr.fetch("1");
    const DiagonalBand band_sel{sel.begin<std::int32_t>(), sel.end<std::int32_t>(), num_bins};

    const auto expected = fetch_pixels<std::int32_t>(sel, num_bins);
    const auto found = band_sel.read_all();

    REQUIRE(expected.size() == 5'288);
    REQUIRE(expected.size() == found.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == found[i]);
    }
  }

  SECTION("gw") {
    constexpr std::uint64_t num_bins = 50;

    const auto sel = clr.fetch();
    const DiagonalBand band_sel{sel.begin<std::int32_t>(), sel.end<std::int32_t>(), num_bins};

    const auto expected = fetch_pixels<std::int32_t>(sel, num_bins);
    const auto found = band_sel.read_all();

    REQUIRE(expected.size() == 54901);
    REQUIRE(expected.size() == found.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == found[i]);
    }
  }

  SECTION("huge num_bins") {
    constexpr auto num_bins = 1'000'000;

    const auto sel = clr.fetch("1");
    const DiagonalBand band_sel{sel.begin<std::int32_t>(), sel.end<std::int32_t>(), num_bins};

    const auto expected = fetch_pixels<std::int32_t>(sel, num_bins);
    const auto found = band_sel.read_all();

    REQUIRE(expected.size() == 5812);
    REQUIRE(expected.size() == found.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == found[i]);
    }
  }

  SECTION("0 num_bins") {
    constexpr std::uint64_t num_bins = 0;

    const auto sel = clr.fetch("1");
    const DiagonalBand band_sel{sel.begin<std::int32_t>(), sel.end<std::int32_t>(), num_bins};

    const auto found = band_sel.read_all();

    REQUIRE(found.empty());
  }

  SECTION("empty range") {
    constexpr std::uint64_t num_bins = 200;

    const auto sel = clr.fetch("1");
    const DiagonalBand band_sel{sel.end<std::int32_t>(), sel.end<std::int32_t>(), num_bins};
    const auto found = band_sel.read_all();

    CHECK(found.empty());
  }
}

TEST_CASE("Transformers (hic): diagonal band", "[transformers][short]") {
  const hic::File hf((datadir / "hic" / "ENCFF993FGR.2500000.hic").string(), 2'500'000);
  SECTION("simple") {
    constexpr std::uint64_t num_bins = 200;

    const auto sel = hf.fetch("chr1");
    const DiagonalBand band_sel{sel.begin<std::int32_t>(), sel.end<std::int32_t>(), num_bins};

    const auto expected = fetch_pixels<std::int32_t>(sel, num_bins);
    const auto found = band_sel.read_all();

    REQUIRE(expected.size() == 4'465);
    REQUIRE(expected.size() == found.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == found[i]);
    }
  }

  SECTION("gw") {
    constexpr std::uint64_t num_bins = 50;

    const auto sel = hf.fetch();
    const DiagonalBand band_sel{sel.begin<std::int32_t>(), sel.end<std::int32_t>(), num_bins};

    const auto expected = fetch_pixels<std::int32_t>(sel, num_bins);
    const auto found = band_sel.read_all();

    REQUIRE(expected.size() == 56'989);
    REQUIRE(expected.size() == found.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == found[i]);
    }
  }

  SECTION("huge num_bins") {
    constexpr auto num_bins = 1'000'000;

    const auto sel = hf.fetch("chr1");
    const DiagonalBand band_sel{sel.begin<std::int32_t>(), sel.end<std::int32_t>(), num_bins};

    const auto expected = fetch_pixels<std::int32_t>(sel, num_bins);
    const auto found = band_sel.read_all();

    REQUIRE(expected.size() == 4'465);
    REQUIRE(expected.size() == found.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(expected[i] == found[i]);
    }
  }

  SECTION("0 num_bins") {
    constexpr std::uint64_t num_bins = 0;

    const auto sel = hf.fetch("chr1");
    const DiagonalBand band_sel{sel.begin<std::int32_t>(), sel.end<std::int32_t>(), num_bins};

    const auto found = band_sel.read_all();

    REQUIRE(found.empty());
  }

  SECTION("empty range") {
    constexpr std::uint64_t num_bins = 200;

    const auto sel = hf.fetch("chr1");
    const DiagonalBand band_sel{sel.end<std::int32_t>(), sel.end<std::int32_t>(), num_bins};
    const auto found = band_sel.read_all();

    CHECK(found.empty());
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::transformers
