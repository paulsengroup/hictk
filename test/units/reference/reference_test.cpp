// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/reference.hpp"

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "hictk/chromosome.hpp"

namespace hictk::test::reference {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Reference", "[reference][short]") {
  // clang-format off
    const std::array<Chromosome, 3> expected{
            Chromosome{0, "chr1", 50001},  //NOLINT
            Chromosome{1, "chr2", 25017},  //NOLINT
            Chromosome{2, "chr3", 10000}   //NOLINT
    };
  // clang-format on

  constexpr std::array<std::string_view, 3> expected_names{"chr1", "chr2", "chr3"};
  constexpr std::array<std::uint32_t, 3> expected_sizes{50001, 25017, 10000};

  SECTION("ctor w/ iterator of chromosomes") {
    const Reference chroms(expected.begin(), expected.end());

    CHECK(chroms.size() == 3);
  }

  SECTION("ctor w/ iterator of chrom names and sizes") {
    const Reference chroms(expected_names.begin(), expected_names.end(), expected_sizes.begin());

    CHECK(chroms.size() == 3);
  }

  SECTION("ctor w/ iterator of chromosomes (duplicates)") {
    std::vector<Chromosome> expected_(expected.begin(), expected.end());
    expected_.push_back(expected.back());

    CHECK_THROWS_WITH(Reference(expected_.begin(), expected_.end()),
                      Catch::Matchers::ContainsSubstring("found multiple entries for chromosome"));
  }

  SECTION("ctor w/ iterator of chrom names and sizes") {
    std::vector<std::string_view> expected_names_(expected_names.begin(), expected_names.end());
    std::vector<std::uint32_t> expected_sizes_(expected_sizes.begin(), expected_sizes.end());

    expected_names_.push_back(expected_names.back());
    expected_sizes_.push_back(expected_sizes.back());

    CHECK_THROWS_WITH(
        Reference(expected_names_.begin(), expected_names_.end(), expected_sizes_.begin()),
        Catch::Matchers::ContainsSubstring("found multiple entries for chromosome"));
  }

  SECTION("contains") {
    const Reference chroms(expected.begin(), expected.end());
    CHECK(chroms.contains(Chromosome{0, "chr1", 50001}));
    CHECK(chroms.contains(0));
    CHECK(chroms.contains("chr1"));

    CHECK_FALSE(chroms.contains(Chromosome{0, "chr0", 50001}));
    CHECK_FALSE(chroms.contains(Chromosome{3, "chr0", 50001}));
    CHECK_FALSE(chroms.contains(7));
    CHECK_FALSE(chroms.contains("chr0"));
    CHECK_FALSE(chroms.contains(""));
  }

  SECTION("at") {
    const Reference chroms(expected.begin(), expected.end());
    CHECK(chroms.at(0) == Chromosome{0, "chr1", 50001});
    CHECK(chroms.at("chr1") == Chromosome{0, "chr1", 50001});

    CHECK_THROWS_AS(chroms.at(3), std::out_of_range);
    CHECK_THROWS_AS(chroms.at("chr0"), std::out_of_range);
  }

  SECTION("opertor[]") {
    const Reference chroms(expected.begin(), expected.end());
    CHECK(chroms[0] == Chromosome{0, "chr1", 50001});
    CHECK(chroms["chr1"] == Chromosome{0, "chr1", 50001});
  }

  SECTION("get_id") {
    const Reference chroms(expected.begin(), expected.end());
    CHECK(chroms.get_id("chr1") == 0);
    CHECK(chroms.get_id("chr3") == 2);

    CHECK_THROWS_AS(chroms.get_id("a"), std::out_of_range);
  }

  SECTION("iteration") {
    const Reference chroms(expected.begin(), expected.end());
    CHECK(std::equal(chroms.begin(), chroms.end(), expected.begin(), expected.end()));
    CHECK(std::equal(chroms.rbegin(), chroms.rend(), expected.rbegin(), expected.rend()));
  }

  SECTION("operators") {
    const Reference chroms1(expected.begin(), expected.end());
    const Reference chroms2(expected.begin(), expected.end() - 1);

    CHECK(chroms1 == chroms1);
    CHECK(chroms1 != chroms2);
  }

  SECTION("accessors") {
    const Reference chroms1(expected.begin(), expected.end());
    const Reference chroms2{{Chromosome{0, "chr1", 1000}, Chromosome{1, "chr123", 5}}};

    CHECK(chroms1.chromosome_with_longest_name().name() == "chr1");
    CHECK(chroms1.longest_chromosome().name() == "chr1");

    CHECK(chroms2.chromosome_with_longest_name().name() == "chr123");
    CHECK(chroms2.longest_chromosome().name() == "chr1");

    const auto& prefix_sum = chroms2.chrom_size_prefix_sum();
    REQUIRE(prefix_sum.size() == 4);
    CHECK(prefix_sum[0] == 0);
    CHECK(prefix_sum[1] == 1000);
    CHECK(prefix_sum[2] == 1005);
    CHECK(prefix_sum[3] == 1006);
  }
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)
}  // namespace hictk::test::reference
