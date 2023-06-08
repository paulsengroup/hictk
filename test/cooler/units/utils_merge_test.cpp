// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "coolerpp/test/self_deleting_folder.hpp"
#include "coolerpp/utils.hpp"

namespace coolerpp::test {
inline const SelfDeletingFolder testdir{true};            // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace coolerpp::test

namespace coolerpp::test::index {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("utils: merge", "[merge][utils][short]") {
  const auto src = datadir / "cooler_test_file.cool";
  const auto dest = testdir() / "cooler_merge_test1.cool";

  const std::array<std::string, 2> sources{src.string(), src.string()};


  SECTION("merge") {
    utils::merge(sources.begin(), sources.end(), dest.string(), true);

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

  SECTION("merge - different resolutions") {
    const auto mclr = datadir / "multires_cooler_test_file.mcool";
    const auto dest1 = testdir() / "cooler_merge_test2.cool";

    const std::array<std::string, 2> sources1{
        fmt::format(FMT_STRING("{}::/resolutions/100000"), mclr.string()),
        fmt::format(FMT_STRING("{}::/resolutions/200000"), mclr.string())};

    CHECK_THROWS_WITH(utils::merge(sources1.begin(), sources1.end(), dest1.string(), true),
                      Catch::Matchers::ContainsSubstring("have different resolutions"));
  }

  SECTION("merge - different reference") {
    const auto clr1 = datadir / "cooler_test_file.cool";
    const auto clr2 = datadir / "ENCFF993FGR.2500000.cool";
    const auto dest2 = testdir() / "cooler_merge_test2.cool";

    const std::array<std::string, 2> sources2{clr1, clr2};

    CHECK_THROWS_WITH(utils::merge(sources2.begin(), sources2.end(), dest2.string(), true),
                      Catch::Matchers::ContainsSubstring("use different reference genomes"));
  }
}
}  // namespace coolerpp::test::index
