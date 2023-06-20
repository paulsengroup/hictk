// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/bin_table.hpp"

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include "hictk/fmt/bin_table.hpp"

namespace hictk::test::bin_table {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
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

    CHECK(fmt::format(FMT_STRING("{}"), bin1) == std::to_string(Bin::null_id));
    CHECK(fmt::format(FMT_STRING("{}"), bin2) == "123");
    CHECK(fmt::format(FMT_STRING("{:bed}"), bin1) == "chr1\t0\t100");
    CHECK(fmt::format(FMT_STRING("{:ucsc}"), bin1) == "chr1:0-100");
    CHECK(fmt::format(FMT_STRING("{:raw}"), bin2) == "123");
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BinTable", "[bin-table][short]") {
  constexpr std::uint32_t bin_size = 5000;
  // clang-format off
  const BinTable table({
      Chromosome{0, "chr1", 50001},
      Chromosome{1, "chr2", 25017},
      Chromosome{2, "chr3", 10000}},
      bin_size);
  // clang-format on

  SECTION("stats") {
    CHECK(BinTable{}.empty());
    CHECK(table.size() == 11 + 6 + 2);
    CHECK(table.num_chromosomes() == 3);
    CHECK(table.bin_size() == bin_size);
  }

  SECTION("at") {
    const auto& chr1 = table.chromosomes().at("chr1");
    const auto& chr2 = table.chromosomes().at("chr2");

    CHECK(table.at(0) == Bin{chr1, 0, bin_size});
    CHECK(table.at(10) == Bin{chr1, 50000, 50001});

    CHECK(table.at(11) == Bin{chr2, 0, bin_size});

    CHECK_THROWS_AS(table.at(table.size()), std::out_of_range);
  }

  SECTION("coord to bin id") {
    const auto& chr2 = table.chromosomes().at("chr2");

    CHECK(table.map_to_bin_id(0, 7500) == 1);
    CHECK(table.map_to_bin_id("chr1", 50000) == 10);
    CHECK(table.map_to_bin_id(chr2, 10) == 11);

    CHECK_THROWS_AS(table.map_to_bin_id("a", 0), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id("chr1", 99999), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id(chr2, 99999), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id(1, 99999), std::out_of_range);
  }

  SECTION("subset") {
    const BinTable expected{{Chromosome{1, "chr2", 25017}}, bin_size};

    CHECK(table.subset(Chromosome{1, "chr2", 25017}) == expected);
    CHECK(table.subset("chr2") == expected);
    CHECK(table.subset(1) == expected);
    CHECK(table.subset("chr1") != expected);

    if constexpr (ndebug_not_defined()) {
      CHECK_THROWS_AS(table.subset(Chromosome{4, "chr5", 1}), std::out_of_range);
    }
    CHECK_THROWS_AS(table.subset("a"), std::out_of_range);
    CHECK_THROWS_AS(table.subset(10), std::out_of_range);
  }

  SECTION("find overlap") {
    const auto& chrom = *table.chromosomes().begin();

    auto its = table.find_overlap({chrom, 10'000, 10'001});
    CHECK(std::distance(its.first, its.second) == 1);

    its = table.find_overlap({chrom, 0, bin_size - 1});
    CHECK(std::distance(its.first, its.second) == 1);

    its = table.find_overlap({chrom, 10'000, 20'000});
    CHECK(std::distance(its.first, its.second) == 2);

    its = table.find_overlap({chrom, 0, chrom.size()});
    const auto table1 = table.subset(chrom);
    CHECK(std::distance(its.first, its.second) == std::distance(table1.begin(), table1.end()));
  }

  SECTION("iterators") {
    const auto& chr1 = table.chromosomes().at("chr1");
    const auto& chr2 = table.chromosomes().at("chr2");
    const auto& chr3 = table.chromosomes().at("chr3");

    // clang-format off
    const std::array<Bin, 19> expected{
       Bin{chr1, 0, 5000},
       Bin{chr1, 5000, 10000},
       Bin{chr1, 10000, 15000},
       Bin{chr1, 15000, 20000},
       Bin{chr1, 20000, 25000},
       Bin{chr1, 25000, 30000},
       Bin{chr1, 30000, 35000},
       Bin{chr1, 35000, 40000},
       Bin{chr1, 40000, 45000},
       Bin{chr1, 45000, 50000},
       Bin{chr1, 50000, 50001},
       Bin{chr2, 0, 5000},
       Bin{chr2, 5000, 10000},
       Bin{chr2, 10000, 15000},
       Bin{chr2, 15000, 20000},
       Bin{chr2, 20000, 25000},
       Bin{chr2, 25000, 25017},
       Bin{chr3, 0, 5000},
       Bin{chr3, 5000, 10000}
    };
    // clang-format on

    REQUIRE(table.size() == expected.size());

    SECTION("forward") {
      auto first_bin = table.begin();
      auto last_bin = table.end();

      // NOLINTNEXTLINE
      for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(*first_bin++ == expected[i]);
      }

      CHECK(first_bin == last_bin);
    }

    SECTION("backward") {
      auto first_bin = table.begin();
      auto last_bin = table.end();

      // NOLINTNEXTLINE
      for (std::size_t i = expected.size(); i != 0; --i) {
        CHECK(*(--last_bin) == expected[i - 1]);
      }

      CHECK(first_bin == last_bin);
    }

    SECTION("operator+") {
      CHECK(table.begin() + 0 == table.begin());
      CHECK(*(table.begin() + 5) == expected[5]);

      auto it = table.begin() + 5;
      for (std::size_t i = 0; i < expected.size() - 5; ++i) {
        CHECK(*(it + i) == expected[i + 5]);  // NOLINT
      }
    }

    SECTION("operator-") {
      CHECK(table.begin() - 0 == table.begin());
      CHECK(*(table.end() - 5) == *(expected.end() - 5));

      auto it1 = table.end();
      auto it2 = expected.end();  // NOLINT
      for (std::size_t i = 1; i < expected.size(); ++i) {
        CHECK(*(it1 - i) == *(it2 - i));  // NOLINT
      }
    }
  }

  SECTION("concretize") {
    const auto concrete_table = table.concretize();

    REQUIRE(concrete_table.chroms.size() == table.size());

    std::size_t i = 0;
    for (const auto& bin : table) {
      CHECK(*concrete_table.chroms[i] == bin.chrom());
      CHECK(concrete_table.bin_starts[i] == bin.start());
      CHECK(concrete_table.bin_ends[i++] == bin.end());
    }
  }
}
}  // namespace hictk::test::bin_table
