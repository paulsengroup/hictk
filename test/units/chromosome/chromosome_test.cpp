// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/chromosome.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <string>
#include <string_view>

#include "hictk/fmt/chromosome.hpp"  // IWYU pragma: keep

namespace hictk::test::chromosome {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Chromosome", "[chromosome][short]") {
  const Chromosome chrom1{0, "chr1", 50001};
  const Chromosome chrom2{1, "chr2", 25017};

  SECTION("accessors") {
    CHECK(chrom1.id() == 0);
    CHECK(chrom1.name() == "chr1");
    CHECK(chrom1.size() == 50001);
    CHECK_FALSE(chrom1.is_all());
    CHECK(Chromosome{0, "All", 10}.is_all());
    CHECK(Chromosome{0, "aLl", 10}.is_all());
    CHECK_FALSE(Chromosome{0, "123", 10}.is_all());
  }

  SECTION("operators") {
    CHECK(!Chromosome{});

    CHECK(chrom1 == chrom1);
    CHECK(chrom1 != chrom2);
    CHECK(chrom1 < chrom2);
    CHECK(chrom1 <= chrom2);
    CHECK(chrom2 > chrom1);
    CHECK(chrom2 >= chrom1);

    CHECK(0U == chrom1);
    CHECK(1U != chrom1);
    CHECK(1U > chrom1);
    CHECK(1U >= chrom1);
    CHECK(0U < chrom2);
    CHECK(0U <= chrom2);

    CHECK("chr1" == chrom1);
    CHECK("chr2" != chrom1);
  }

  SECTION("fmt") {
    CHECK(fmt::format(FMT_STRING("{}"), chrom1) == "chr1:50001");
    CHECK(fmt::format(FMT_STRING("{:tsv}"), chrom1) == "chr1\t50001");
    CHECK(fmt::format(FMT_STRING("{:ucsc}"), chrom1) == "chr1:50001");
  }
}

}  // namespace hictk::test::chromosome
