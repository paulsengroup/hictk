// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/bin_table.hpp"  // IWYU pragma: keep

#include <fmt/format.h>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>

#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/genomic_interval.hpp"

namespace hictk::test::bin_table {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BinTable (fixed bins)", "[bin-table][short]") {
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
    CHECK(table.resolution() == bin_size);
  }

  SECTION("at") {
    const auto& chr1 = table.chromosomes().at("chr1");
    const auto& chr2 = table.chromosomes().at("chr2");

    CHECK(table.at(0) == Bin{chr1, 0, bin_size});
    CHECK(table.at(10) == Bin{chr1, 50000, 50001});
    CHECK(table.at(11) == Bin{chr2, 0, bin_size});

    CHECK(table.at(chr1, bin_size - 1).id() == 0);
    CHECK(table.at(chr1, 50000).id() == 10);
    CHECK(table.at(chr2, 1).id() == 11);

    CHECK_THROWS_AS(table.at(table.size()), std::out_of_range);
    CHECK_THROWS_AS(table.at(chr1, 50001), std::out_of_range);
    CHECK_THROWS_AS(table.at(chr2, 26000), std::out_of_range);
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
    const BinTable expected{Reference{Chromosome{1, "chr2", 25017}}, bin_size};

    CHECK(table.subset(Chromosome{1, "chr2", 25017}) == expected);
    CHECK(table.subset("chr2") == expected);
    CHECK(table.subset(1) == expected);
    CHECK(table.subset("chr1") != expected);
    CHECK(table.subset("chr2").subset("chr2") == expected);

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

  SECTION("accessors") {
    CHECK(table.type() == BinTable::Type::fixed);
    CHECK_NOTHROW(table.get<BinTableFixed>());
    CHECK_THROWS(table.get<BinTableVariable<>>());
  }

  SECTION("operator==") {
    CHECK(BinTable(table.chromosomes(), 10) == BinTable(table.chromosomes(), 10));
    CHECK(BinTable(table.chromosomes(), 10) != BinTable(table.chromosomes(), 20));
    CHECK(BinTable(Reference{table.chromosomes().begin(), table.chromosomes().end() - 1}, 10) !=
          BinTable(table.chromosomes(), 10));
  }

  SECTION("iterators") {
    const auto& chr1 = table.chromosomes().at("chr1");
    const auto& chr2 = table.chromosomes().at("chr2");
    const auto& chr3 = table.chromosomes().at("chr3");

    // clang-format off
    const std::array<Bin, 19> expected{
       Bin{0,   0, chr1, 0, 5000},
       Bin{1,   1, chr1, 5000, 10000},
       Bin{2,   2, chr1, 10000, 15000},
       Bin{3,   3, chr1, 15000, 20000},
       Bin{4,   4, chr1, 20000, 25000},
       Bin{5,   5, chr1, 25000, 30000},
       Bin{6,   6, chr1, 30000, 35000},
       Bin{7,   7, chr1, 35000, 40000},
       Bin{8,   8, chr1, 40000, 45000},
       Bin{9,   9, chr1, 45000, 50000},
       Bin{10, 10, chr1, 50000, 50001},
       Bin{11,  0, chr2, 0, 5000},
       Bin{12,  1, chr2, 5000, 10000},
       Bin{13,  2, chr2, 10000, 15000},
       Bin{14,  3, chr2, 15000, 20000},
       Bin{15,  4, chr2, 20000, 25000},
       Bin{16,  5, chr2, 25000, 25017},
       Bin{17,  0, chr3, 0, 5000},
       Bin{18,  1, chr3, 5000, 10000}
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

      CHECK_THROWS_AS(first_bin++, std::out_of_range);
      CHECK_THROWS_AS(last_bin++, std::out_of_range);
    }

    SECTION("backward") {
      auto first_bin = table.begin();
      auto last_bin = table.end();

      // NOLINTNEXTLINE
      for (std::size_t i = expected.size(); i != 0; --i) {
        CHECK(*(--last_bin) == expected[i - 1]);
      }

      CHECK(first_bin == last_bin);

      CHECK_THROWS_AS(--first_bin, std::out_of_range);
      CHECK_THROWS_AS(--last_bin, std::out_of_range);
    }

    SECTION("operator+") {
      CHECK(table.begin() + 0 == table.begin());
      CHECK(*(table.begin() + 5) == expected[5]);

      auto it = table.begin() + 5;
      for (std::size_t i = 0; i < expected.size() - 5; ++i) {
        CHECK(*(it + i) == expected[i + 5]);  // NOLINT
      }

      CHECK_THROWS_AS(it + 100, std::out_of_range);
    }

    SECTION("operator-") {
      CHECK(table.begin() - 0 == table.begin());
      CHECK(*(table.end() - 5) == *(expected.end() - 5));

      auto it1 = table.end();
      auto it2 = expected.end();  // NOLINT
      for (std::size_t i = 1; i < expected.size(); ++i) {
        CHECK(*(it1 - i) == *(it2 - i));  // NOLINT
      }

      CHECK_THROWS_AS(it1 - 100, std::out_of_range);
    }

    SECTION("accessors") {
      CHECK_NOTHROW(table.begin().get<BinTableFixed::iterator>());
      CHECK_THROWS(table.begin().get<BinTableVariable<>::iterator>());
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("BinTable (variable bins)", "[bin-table][short]") {
  const Chromosome chrom1{0, "chr1", 32};
  const Chromosome chrom2{1, "chr2", 32};

  const std::vector<std::uint32_t> start_pos{0, 8, 15, 23, 0, 5, 10, 26};
  const std::vector<std::uint32_t> end_pos{8, 15, 23, 32, 5, 10, 26, 32};

  const BinTable table({chrom1, chrom2}, start_pos, end_pos);

  SECTION("stats") {
    CHECK(BinTable{}.empty());
    CHECK(table.size() == start_pos.size());
    CHECK(table.num_chromosomes() == 2);
    CHECK(table.resolution() == 0);
  }

  SECTION("at") {
    const auto& chr1 = table.chromosomes().at("chr1");
    const auto& chr2 = table.chromosomes().at("chr2");

    CHECK(table.at(0) == Bin{chr1, 0, 8});
    CHECK(table.at(3) == Bin{chr1, 23, 32});
    CHECK(table.at(4) == Bin{chr2, 0, 5});

    CHECK(table.at(chr1, 0).id() == 0);
    CHECK(table.at(chr1, 7).id() == 0);
    CHECK(table.at(chr1, 8).id() == 1);

    CHECK(table.at(chr1, 23).id() == 3);
    CHECK(table.at(chr2, 4).id() == 4);

    CHECK_THROWS_AS(table.at(table.size()), std::out_of_range);
    CHECK_THROWS_AS(table.at(chr1, 32), std::out_of_range);
    CHECK_THROWS_AS(table.at(chr2, 32), std::out_of_range);
  }

  SECTION("coord to bin id") {
    const auto& chr2 = table.chromosomes().at("chr2");

    CHECK(table.map_to_bin_id(0, 8) == 1);
    CHECK(table.map_to_bin_id("chr1", 25) == 3);
    CHECK(table.map_to_bin_id(chr2, 9) == 5);

    CHECK_THROWS_AS(table.map_to_bin_id("a", 0), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id("chr1", 33), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id(chr2, 50), std::out_of_range);
    CHECK_THROWS_AS(table.map_to_bin_id(1, 50), std::out_of_range);
  }

  SECTION("subset") {
    const std::vector<std::uint32_t> start_pos_{0, 5, 10, 26};
    const std::vector<std::uint32_t> end_pos_{5, 10, 26, 32};
    const BinTable expected{Reference{Chromosome{1, "chr2", 32}}, start_pos_, end_pos_};

    CHECK(table.subset(Chromosome{1, "chr2", 32}) == expected);
    CHECK(table.subset("chr2") == expected);
    CHECK(table.subset(1) == expected);
    CHECK(table.subset("chr1") != expected);
    CHECK(table.subset("chr2").subset("chr2") == expected);

    if constexpr (ndebug_not_defined()) {
      CHECK_THROWS_AS(table.subset(Chromosome{4, "chr5", 1}), std::out_of_range);
    }
    CHECK_THROWS_AS(table.subset("a"), std::out_of_range);
    CHECK_THROWS_AS(table.subset(10), std::out_of_range);
  }

  SECTION("find overlap") {
    const auto& chrom = *table.chromosomes().begin();

    auto its = table.find_overlap({chrom, 8, 9});
    CHECK(std::distance(its.first, its.second) == 1);

    its = table.find_overlap({chrom, 8, 15 - 1});
    CHECK(std::distance(its.first, its.second) == 1);

    its = table.find_overlap({chrom, 14, 23});
    CHECK(std::distance(its.first, its.second) == 2);

    its = table.find_overlap({chrom, 0, chrom.size()});
    CHECK(std::distance(its.first, its.second) == 4);
  }

  SECTION("accessors") {
    CHECK(table.type() == BinTable::Type::variable);
    CHECK_NOTHROW(table.get<BinTableVariable<>>());
    CHECK_THROWS(table.get<BinTableFixed>());
  }

  SECTION("invalid bins") {
    SECTION("bins out of order") {
      const std::vector<std::uint32_t> start_pos1{0, 8, 7};
      const std::vector<std::uint32_t> end_pos1{8, 15, 23};

      CHECK_THROWS_WITH(BinTable({chrom1, chrom2}, start_pos1, end_pos1),
                        Catch::Matchers::ContainsSubstring("not sorted"));

      const std::vector<std::uint32_t> start_pos2{0, 8, 15};
      const std::vector<std::uint32_t> end_pos2{8, 15, 14};

      CHECK_THROWS_WITH(BinTable({chrom1, chrom2}, start_pos2, end_pos2),
                        Catch::Matchers::ContainsSubstring("not sorted"));
    }
    SECTION("gap between bins") {
      const std::vector<std::uint32_t> start_pos1{0, 8, 16};
      const std::vector<std::uint32_t> end_pos1{8, 15, 23};

      CHECK_THROWS_WITH(BinTable({chrom1, chrom2}, start_pos1, end_pos1),
                        Catch::Matchers::ContainsSubstring("gap between bins"));

      const std::vector<std::uint32_t> start_pos2{1, 8, 16};
      const std::vector<std::uint32_t> end_pos2{8, 15, 23};

      CHECK_THROWS_WITH(BinTable({chrom1, chrom2}, start_pos2, end_pos2),
                        Catch::Matchers::ContainsSubstring("does not start from zero"));
    }

    SECTION("start pos >= end pos") {
      const std::vector<std::uint32_t> start_pos1{0, 8, 10, 15};
      const std::vector<std::uint32_t> end_pos1{0, 10, 15, 23};

      CHECK_THROWS_WITH(BinTable({chrom1, chrom2}, start_pos1, end_pos1),
                        Catch::Matchers::ContainsSubstring("start_pos >= end_pos"));
    }

    SECTION("number of chromosome mismatch") {
      const Chromosome chrom3{2, "chr3", 32};
      CHECK_THROWS_WITH(BinTable({chrom1, chrom2, chrom3}, start_pos, end_pos),
                        Catch::Matchers::ContainsSubstring("unexpected number of chromosomes"));
    }
  }

  SECTION("operator==") {
    CHECK(BinTable(table.chromosomes(), start_pos, end_pos) ==
          BinTable(table.chromosomes(), start_pos, end_pos));

    const std::vector<std::uint32_t> start_pos1{0, 0};
    const std::vector<std::uint32_t> end_pos1{32, 32};
    CHECK(BinTable(table.chromosomes(), start_pos, end_pos) !=
          BinTable(table.chromosomes(), start_pos1, end_pos1));

    const std::vector<std::uint32_t> start_pos2{0};
    const std::vector<std::uint32_t> end_pos2{32};
    CHECK(BinTable(Reference{table.chromosomes().begin(), table.chromosomes().end() - 1},
                   start_pos2, end_pos2) != BinTable(table.chromosomes(), start_pos, end_pos));

    CHECK(BinTable(table.chromosomes(), start_pos, end_pos) != BinTable(table.chromosomes(), 10));
  }

  SECTION("iterators") {
    const auto& chr1 = table.chromosomes().at("chr1");
    const auto& chr2 = table.chromosomes().at("chr2");

    // clang-format off
    const std::array<Bin, 8> expected{
       Bin{0, 0, chr1, 0,  8},
       Bin{1, 1, chr1, 8,  15},
       Bin{2, 2, chr1, 15, 23},
       Bin{3, 3, chr1, 23, 32},
       Bin{4, 0, chr2, 0,  5},
       Bin{5, 1, chr2, 5,  10},
       Bin{6, 2, chr2, 10, 26},
       Bin{7, 3, chr2, 26, 32}
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

      CHECK_THROWS_AS(first_bin++, std::out_of_range);
      CHECK_THROWS_AS(last_bin++, std::out_of_range);
    }

    SECTION("backward") {
      auto first_bin = table.begin();
      auto last_bin = table.end();

      // NOLINTNEXTLINE
      for (std::size_t i = expected.size(); i != 0; --i) {
        CHECK(*(--last_bin) == expected[i - 1]);
      }

      CHECK(first_bin == last_bin);

      CHECK_THROWS_AS(--first_bin, std::out_of_range);
      CHECK_THROWS_AS(--last_bin, std::out_of_range);
    }

    SECTION("operator+") {
      CHECK(table.begin() + 0 == table.begin());
      CHECK(*(table.begin() + 5) == expected[5]);

      auto it = table.begin() + 5;
      for (std::size_t i = 0; i < expected.size() - 5; ++i) {
        CHECK(*(it + i) == expected[i + 5]);  // NOLINT
      }

      CHECK_THROWS_AS(it + 100, std::out_of_range);
    }

    SECTION("operator-") {
      CHECK(table.begin() - 0 == table.begin());
      CHECK(*(table.end() - 5) == *(expected.end() - 5));

      auto it1 = table.end();
      auto it2 = expected.end();  // NOLINT
      for (std::size_t i = 1; i < expected.size(); ++i) {
        CHECK(*(it1 - i) == *(it2 - i));  // NOLINT
      }

      CHECK_THROWS_AS(it1 - 100, std::out_of_range);
    }

    SECTION("accessors") {
      CHECK_NOTHROW(table.begin().get<BinTableVariable<>::iterator>());
      CHECK_THROWS(table.begin().get<BinTableFixed::iterator>());
    }
  }
}

}  // namespace hictk::test::bin_table
