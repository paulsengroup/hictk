// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <string>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::utils {
static const auto& testdir = hictk::test::testdir;
static const auto& datadir = hictk::test::datadir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: utils merge", "[merge][utils][long]") {
  SECTION("merge int") {
    const auto src = datadir / "cooler" / "cooler_test_file.cool";
    const auto dest = testdir() / "cooler_merge_test_int.cool";

    const std::array<std::string, 2> sources{src.string(), src.string()};

    cooler::utils::merge<std::int32_t>(sources.begin(), sources.end(), dest.string(), true, 1'000);

    const auto clr1 = File::open_read_once(src.string());
    const auto clr2 = File::open_read_once(dest.string());

    auto first1 = clr1.begin<std::int32_t>();
    auto last1 = clr1.end<std::int32_t>();

    auto first2 = clr2.begin<std::int32_t>();
    auto last2 = clr2.end<std::int32_t>();

    REQUIRE(std::distance(first1, last1) == std::distance(first2, last2));
    while (first1 != last1) {
      CHECK(first1->bin1_id == first2->bin1_id);
      CHECK(first1->bin2_id == first2->bin2_id);
      CHECK(2 * first1->count == first2->count);
      ++first1;
      ++first2;
    }
  }

  SECTION("merge float") {
    const auto src = datadir / "cooler" / "cooler_test_file_float.cool";
    const auto dest = testdir() / "cooler_merge_test_float.cool";

    const std::array<std::string, 2> sources{src.string(), src.string()};

    cooler::utils::merge<double>(sources.begin(), sources.end(), dest.string(), true, 1'000);

    const auto clr1 = File::open_read_once(src.string());
    const auto clr2 = File::open_read_once(dest.string());

    auto first1 = clr1.begin<std::int32_t>();
    auto last1 = clr1.end<std::int32_t>();

    auto first2 = clr2.begin<std::int32_t>();
    auto last2 = clr2.end<std::int32_t>();

    REQUIRE(std::distance(first1, last1) == std::distance(first2, last2));
    while (first1 != last1) {
      CHECK(first1->bin1_id == first2->bin1_id);
      CHECK(first1->bin2_id == first2->bin2_id);
      CHECK(2 * first1->count == first2->count);
      ++first1;
      ++first2;
    }
  }

  SECTION("merge chromosomes") {
    const auto src = datadir / "cooler" / "cooler_test_file.cool";
    const auto dest = testdir() / "cooler_merge_test2.cool";
    std::vector<std::string> sources{};
    {
      const cooler::File clr(src.string());

      for (const auto& chrom : clr.chromosomes()) {
        sources.emplace_back((testdir() / std::string{chrom.name()}).string());

        auto clr1 = cooler::File::create(sources.back(), clr.chromosomes(), clr.resolution());
        const auto sel = clr.fetch(chrom.name());
        clr1.append_pixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>());
      }
    }

    cooler::utils::merge<std::int32_t>(sources.begin(), sources.end(), dest.string(), true, 1000);

    const auto clr1 = File::open_read_once(src.string());
    const auto clr2 = File::open_read_once(dest.string());

    for (const auto& chrom : clr1.chromosomes()) {
      auto sel1 = clr1.fetch(chrom.name());
      auto sel2 = clr2.fetch(chrom.name());

      auto first1 = sel1.begin<std::int32_t>();
      auto last1 = sel1.end<std::int32_t>();
      auto first2 = sel2.begin<std::int32_t>();
      auto last2 = sel2.end<std::int32_t>();

      REQUIRE(std::distance(first1, last1) == std::distance(first2, last2));
      while (first1 != last1) {
        CHECK(first1->bin1_id == first2->bin1_id);
        CHECK(first1->bin2_id == first2->bin2_id);
        CHECK(first1->count == first2->count);
        ++first1;
        ++first2;
      }
    }
  }

  SECTION("merge - different resolutions") {
    const auto mclr = datadir / "cooler" / "multires_cooler_test_file.mcool";
    const auto dest1 = testdir() / "cooler_merge_test3.cool";

    const std::array<std::string, 2> sources1{
        fmt::format(FMT_STRING("{}::/resolutions/100000"), mclr.string()),
        fmt::format(FMT_STRING("{}::/resolutions/200000"), mclr.string())};

    CHECK_THROWS_WITH(
        cooler::utils::merge<std::int32_t>(sources1.begin(), sources1.end(), dest1.string(), true),
        Catch::Matchers::ContainsSubstring("have different resolutions"));
  }

  SECTION("merge - different reference") {
    const auto clr1 = datadir / "cooler" / "cooler_test_file.cool";
    const auto clr2 = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const auto dest2 = testdir() / "cooler_merge_test2.cool";

    const std::array<std::string, 2> sources2{clr1.string(), clr2.string()};

    CHECK_THROWS_WITH(
        cooler::utils::merge<std::int32_t>(sources2.begin(), sources2.end(), dest2.string(), true),
        Catch::Matchers::ContainsSubstring("use different reference genomes"));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::utils
