// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "hictk/cooler/utils.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk::test {
inline const internal::TmpDir testdir{true};                     // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::cooler::test::utils {
inline const auto& testdir = hictk::test::testdir;
inline const auto& datadir = hictk::test::datadir;

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: utils merge", "[merge][utils][long]") {
  const auto src = datadir / "cooler_test_file.cool";

  SECTION("merge with self") {
    const auto dest = testdir() / "cooler_merge_test1.cool";
    const std::array<std::string, 2> sources{src.string(), src.string()};
    cooler::utils::merge(sources.begin(), sources.end(), dest.string(), true, 1000);

    const auto clr1 = File::open_read_only_read_once(src.string());
    const auto clr2 = File::open_read_only_read_once(dest.string());

    auto first1 = clr1.begin<std::int32_t>();
    auto last1 = clr1.end<std::int32_t>();

    auto first2 = clr2.begin<std::int32_t>();
    auto last2 = clr2.end<std::int32_t>();

    REQUIRE(std::distance(first1, last1) == std::distance(first2, last2));
    while (first1 != last1) {
      CHECK(first1->coords == first2->coords);
      CHECK(2 * first1->count == first2->count);
      ++first1;
      ++first2;
    }
  }

  SECTION("merge chromosomes") {
    const auto dest = testdir() / "cooler_merge_test2.cool";
    std::vector<std::string> sources{};
    {
      auto clr = cooler::File::open_read_only(src.string());

      for (const auto& chrom : clr.chromosomes()) {
        sources.emplace_back((testdir() / std::string{chrom.name()}).string());

        auto clr1 =
            cooler::File::create_new_cooler(sources.back(), clr.chromosomes(), clr.bin_size());
        const auto sel = clr.fetch<std::int32_t>(chrom.name());
        clr1.append_pixels(sel.begin(), sel.end());
      }
    }

    cooler::utils::merge(sources.begin(), sources.end(), dest.string(), true, 1000);

    const auto clr1 = File::open_read_only_read_once(src.string());
    const auto clr2 = File::open_read_only_read_once(dest.string());

    for (const auto& chrom : clr1.chromosomes()) {
      auto sel1 = clr1.fetch<std::int32_t>(chrom.name());
      auto sel2 = clr2.fetch<std::int32_t>(chrom.name());

      auto first1 = sel1.begin();
      auto last1 = sel1.end();
      auto first2 = sel2.begin();
      auto last2 = sel2.end();

      REQUIRE(std::distance(first1, last1) == std::distance(first2, last2));
      while (first1 != last1) {
        CHECK(first1->coords == first2->coords);
        CHECK(first1->count == first2->count);
        ++first1;
        ++first2;
      }
    }
  }

  SECTION("merge - different resolutions") {
    const auto mclr = datadir / "multires_cooler_test_file.mcool";
    const auto dest1 = testdir() / "cooler_merge_test3.cool";

    const std::array<std::string, 2> sources1{
        fmt::format(FMT_STRING("{}::/resolutions/100000"), mclr.string()),
        fmt::format(FMT_STRING("{}::/resolutions/200000"), mclr.string())};

    CHECK_THROWS_WITH(cooler::utils::merge(sources1.begin(), sources1.end(), dest1.string(), true),
                      Catch::Matchers::ContainsSubstring("have different resolutions"));
  }

  SECTION("merge - different reference") {
    const auto clr1 = datadir / "cooler_test_file.cool";
    const auto clr2 = datadir / "ENCFF993FGR.2500000.cool";
    const auto dest2 = testdir() / "cooler_merge_test2.cool";

    const std::array<std::string, 2> sources2{clr1.string(), clr2.string()};

    CHECK_THROWS_WITH(cooler::utils::merge(sources2.begin(), sources2.end(), dest2.string(), true),
                      Catch::Matchers::ContainsSubstring("use different reference genomes"));
  }
}
}  // namespace hictk::cooler::test::utils
