// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/index.hpp"

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"

namespace hictk::cooler::test::index {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: index ctor", "[index][short]") {
  constexpr std::uint32_t bin_size = 100;
  const auto bins = std::make_shared<const BinTable>(
      Reference{Chromosome{0, "chr1", 10001}, Chromosome{1, "chr2", 5000}}, bin_size);

  const Index idx(bins);

  CHECK(idx.resolution() == bin_size);
  CHECK(idx.chromosomes().size() == 2);
  CHECK(idx.size() == 151);

  CHECK(idx.size("chr1") == 101);
  CHECK(idx.size(0) == 101);

  CHECK(idx.size("chr2") == 50);
  CHECK(idx.size(1) == 50);

  CHECK_THROWS(idx.size("chr3"));
  CHECK_THROWS(idx.size(99));
}

TEST_CASE("Cooler: index accessors", "[index][short]") {
  constexpr std::uint32_t bin_size = 100;
  const auto bins = std::make_shared<const BinTable>(
      Reference{Chromosome{0, "chr1", 10001}, Chromosome{1, "chr2", 5000}}, bin_size);

  const Index idx(bins);

  CHECK(Index{}.empty());
  CHECK(Index{}.bins().empty());
  CHECK(Index{}.chromosomes().empty());
  CHECK(Index{}.resolution() == 0);

  CHECK(idx.empty("chr1"));
  CHECK_NOTHROW(idx.at("chr1"));
  CHECK_THROWS_AS(idx.at("chr123"), std::out_of_range);

  CHECK(idx.contains("chr1"));
  CHECK_FALSE(idx.contains("chr123"));
}

TEST_CASE("Cooler: index offset setters and getters", "[index][short]") {
  constexpr std::uint32_t bin_size = 10;
  const Chromosome chrom1{0, "chr1", 100};
  const auto bins = std::make_shared<const BinTable>(Reference{chrom1}, bin_size);

  constexpr auto fill_value = std::numeric_limits<std::uint64_t>::max();

  Index idx(bins);

  SECTION("by pos") {
    idx.set_offset_by_pos(chrom1, 22, 1);
    idx.set_offset_by_pos(0, 55, 1);

    for (std::uint32_t pos = 0; pos < 100; ++pos) {
      const auto row_idx = conditional_static_cast<std::size_t>(pos / bin_size);
      const std::uint64_t expected = (row_idx == 2 || row_idx == 5) ? 1 : fill_value;

      CHECK(idx.get_offset_by_row_idx(0, row_idx) == expected);

      CHECK(idx.get_offset_by_pos(chrom1, pos) == expected);
      CHECK(idx.get_offset_by_pos(0, pos) == expected);
    }
  }

  SECTION("by row idx") {
    idx.set_offset_by_row_idx(0, 2, 1);
    idx.set_offset_by_row_idx(0, 5, 1);

    for (std::uint32_t pos = 0; pos < 100; ++pos) {
      const auto row_idx = conditional_static_cast<std::size_t>(pos / bin_size);
      const std::size_t expected = (row_idx == 2 || row_idx == 5) ? 1 : fill_value;

      CHECK(idx.get_offset_by_row_idx(0, row_idx) == expected);

      CHECK(idx.get_offset_by_pos(chrom1, pos) == expected);
      CHECK(idx.get_offset_by_pos(0, pos) == expected);
    }
  }

  SECTION("by bin ID") {
    idx.set_offset_by_bin_id(9, 9);
    CHECK(idx.get_offset_by_pos(chrom1, 99) == 9);
    CHECK(idx.get_offset_by_bin_id(9) == 9);
  }

  SECTION("out of bound access") {
    CHECK_THROWS_WITH(idx.get_offset_by_pos(chrom1, 999),
                      Catch::Matchers::ContainsSubstring("row maps outside of chromosome"));

    CHECK_THROWS_WITH(idx.get_offset_by_row_idx(0, 999),
                      Catch::Matchers::ContainsSubstring("row maps outside of chromosome"));
  }
}

TEST_CASE("Cooler: index iterator", "[index][short]") {
  constexpr std::uint32_t bin_size = 1000;
  const auto bins = std::make_shared<const BinTable>(
      Reference{Chromosome{0, "chr1", 10001}, Chromosome{1, "chr2", 5000}}, bin_size);

  // Assume there are 10 pixels per row
  constexpr std::array<std::size_t, 11> chr1_offsets{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  constexpr std::array<std::size_t, 5> chr2_offsets{110, 120, 130, 140, 150};

  Index idx(bins);
  for (std::size_t i = 0; i < chr1_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(0, i, chr1_offsets[i]);
  }

  for (std::size_t i = 0; i < chr2_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(1, i, chr2_offsets[i]);
  }

  auto first_offset = idx.begin();
  const auto last_offset = idx.end();
  REQUIRE(first_offset != last_offset);

  for (std::size_t i = 0; i < idx.size(); ++i) {
    CHECK(*first_offset++ == conditional_static_cast<std::uint64_t>(i * 10));
  }

  REQUIRE(first_offset != last_offset);
  CHECK(*first_offset == 0);
  idx.finalize(16);
  CHECK(*first_offset++ == idx.size());
  CHECK(first_offset == last_offset);
}

TEST_CASE("Cooler: index validation", "[index][short]") {
  constexpr std::uint32_t bin_size = 1000;
  const auto bins = std::make_shared<const BinTable>(
      Reference{Chromosome{0, "chr1", 10001}, Chromosome{1, "chr2", 5000}}, bin_size);

  // Assume there are 10 pixels per row
  constexpr std::array<std::size_t, 11> chr1_offsets{0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  constexpr std::array<std::size_t, 5> chr2_offsets{110, 120, 130, 140, 150};

  Index idx(bins);
  for (std::size_t i = 0; i < chr1_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(0, i, chr1_offsets[i]);
  }

  for (std::size_t i = 0; i < chr2_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(1, i, chr2_offsets[i]);
  }

  SECTION("valid index") { CHECK_NOTHROW(idx.validate()); }
  SECTION("first offset is not zero") {
    idx.set_offset_by_row_idx(0, 0, 1);
    CHECK_THROWS_WITH(idx.validate(),
                      Catch::Matchers::ContainsSubstring("first offset is not zero"));
  }
  SECTION("offsets for adjacent chromosomes are not in ascending order") {
    idx.set_offset_by_row_idx(1, 0, 99);
    CHECK_THROWS_WITH(idx.validate(), Catch::Matchers::ContainsSubstring(
                                          "offset for bin chr2:0-1000 should be >= 100, found 99"));
  }

  SECTION("offsets are not sorted") {
    idx.set_offset_by_row_idx(1, 2, 150);
    CHECK_THROWS_WITH(idx.validate(),
                      Catch::Matchers::ContainsSubstring("offsets are not in ascending order"));
  }
}

TEST_CASE("Cooler: index compute chromosome offsets", "[index][short]") {
  constexpr std::uint32_t bin_size = 1000;
  const auto bins = std::make_shared<const BinTable>(
      Reference{Chromosome{0, "chr1", 10001}, Chromosome{1, "chr2", 5000}}, bin_size);

  // Assume there are 10 pixels per row
  constexpr std::array<std::size_t, 11> chr1_offsets{10, 20, 30, 40, 50, 60, 70, 80, 90, 100};
  constexpr std::array<std::size_t, 5> chr2_offsets{110, 120, 130, 140, 150};

  Index idx(bins);

  for (std::size_t i = 0; i < chr1_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(0, i, chr1_offsets[i]);
  }

  for (std::size_t i = 0; i < chr2_offsets.size(); ++i) {
    idx.set_offset_by_row_idx(1, i, chr2_offsets[i]);
  }

  const auto chrom_offsets = idx.compute_chrom_offsets();
  REQUIRE(chrom_offsets.size() == bins->num_chromosomes() + 1);

  CHECK(chrom_offsets[0] == 0);
  CHECK(chrom_offsets[1] == chr1_offsets.size());
  CHECK(chrom_offsets[2] == chr1_offsets.size() + chr2_offsets.size());
}
// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::index
