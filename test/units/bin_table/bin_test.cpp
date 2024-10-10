// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/bin.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <stdexcept>

#include "hictk/fmt/bin.hpp"

namespace hictk::test::bin_table {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Bin", "[bin][short]") {
  const Chromosome chrom1{0, "chr1", 50};
  const Chromosome chrom2{1, "chr2", 10};
  SECTION("Ctors") {
    CHECK(Bin{chrom1, 1, 2}.has_null_id());
    CHECK_FALSE(Bin{0, 0, chrom1, 1, 2}.has_null_id());
    CHECK_FALSE(Bin{0, 0, GenomicInterval{chrom1, 1, 2}}.has_null_id());
  }

  SECTION("Accessors") {
    const Bin bin1{chrom1, 1, 2};
    const Bin bin2{10, 5, chrom1, 1, 2};

    CHECK(bin1.id() == Bin::null_id);
    CHECK(bin2.id() == 10);
    CHECK(bin2.rel_id() == 5);

    CHECK(bin1.interval() == GenomicInterval{chrom1, 1, 2});

    CHECK(bin2.chrom() == chrom1);
    CHECK(bin2.start() == 1);
    CHECK(bin2.end() == 2);
  }

  SECTION("operators (wo/ id)") {
    const Bin bin0{};
    const Bin bin1{chrom1, 1, 2};
    const Bin bin2{chrom1, 2, 3};

    const Bin bin3{chrom2, 1, 2};

    CHECK(!bin0);
    CHECK(!!bin1);

    CHECK(bin1 != bin2);
    CHECK(bin1 != bin3);

    CHECK(bin1 < bin2);
    CHECK(bin1 < bin3);

    CHECK(bin1 <= bin2);
    CHECK(bin1 <= bin3);

    CHECK(bin2 > bin1);
    CHECK(bin3 > bin1);

    CHECK(bin2 >= bin1);
    CHECK(bin3 >= bin1);
  }

  SECTION("operators (w/ id)") {
    const Bin bin1{0, 0, chrom1, 1, 2};
    const Bin bin2{1, 1, chrom1, 2, 3};

    const Bin bin3{10, 10, chrom2, 1, 2};
    const Bin bin4{10, 10, chrom2, 10, 20};

    CHECK(bin1 != bin2);
    CHECK(bin1 != bin3);

    // This is true because they have the same ID.
    // However, comparing bins with same ID and different interval is UB
    CHECK(bin3 == bin4);

    CHECK(bin1 < bin2);
    CHECK(bin1 < bin3);

    CHECK(bin1 <= bin2);
    CHECK(bin1 <= bin3);

    CHECK(bin2 > bin1);
    CHECK(bin3 > bin1);

    CHECK(bin2 >= bin1);
    CHECK(bin3 >= bin1);
  }

  SECTION("fmt") {
    const Bin bin1{chrom1, 0, 100};
    const Bin bin2{123, 123, chrom1, 0, 100};

    CHECK(fmt::format(FMT_STRING("{}"), bin1) == fmt::to_string(Bin::null_id));
    CHECK(fmt::format(FMT_STRING("{}"), bin2) == "123");
    CHECK(fmt::format(FMT_STRING("{:bed}"), bin1) == "chr1\t0\t100");
    CHECK(fmt::format(FMT_STRING("{:ucsc}"), bin1) == "chr1:0-100");
    CHECK(fmt::format(FMT_STRING("{:raw}"), bin2) == "123");
  }
}
// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::bin_table
