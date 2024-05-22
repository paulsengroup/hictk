// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <numeric>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::pixel_selector {

template <typename N>
static std::ptrdiff_t generate_test_data(const std::filesystem::path& path, const Reference& chroms,
                                         std::uint32_t bin_size) {
  auto f = File::create<N>(path.string(), chroms, bin_size, true);

  const auto num_bins = f.bins().size();

  std::vector<ThinPixel<N>> pixels;

  N n = 1;
  for (std::uint64_t bin1_id = 0; bin1_id < num_bins; ++bin1_id) {
    for (std::uint64_t bin2_id = bin1_id; bin2_id < num_bins; ++bin2_id) {
      pixels.emplace_back(ThinPixel<N>{bin1_id, bin2_id, n++});
    }
  }
  f.append_pixels(pixels.begin(), pixels.end());
  return static_cast<std::ptrdiff_t>(pixels.size());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: pixel selector 1D queries", "[pixel_selector][short]") {
  const auto path1 = testdir() / "pixel_selector_devel.cool";

  const Reference chroms{Chromosome{0, "chr1", 1000}, Chromosome{1, "chr2", 100}};
  constexpr std::uint32_t bin_size = 10;
  using T = std::uint32_t;

  const auto expected_nnz = generate_test_data<T>(path1, chroms, bin_size);

  const File f(path1.string());
  REQUIRE(std::distance(f.begin<T>(), f.end<T>()) == expected_nnz);

  SECTION("query overlaps chrom start") {
    auto selector = f.fetch("chr1:0-20");
    const auto pixels = selector.read_all<T>();
    REQUIRE(pixels.size() == 3);

    CHECK(pixels[0].count == 1);
    CHECK(pixels[1].count == 2);
    CHECK(pixels[2].count == 111);
  }

  SECTION("query overlaps chrom end") {
    auto selector = f.fetch("chr1:980-1000");
    const auto pixels = selector.read_all<T>();
    REQUIRE(pixels.size() == 3);

    CHECK(pixels[0].count == 6028);
    CHECK(pixels[1].count == 6029);
    CHECK(pixels[2].count == 6040);
  }

  SECTION("query does not overlap chrom boundaries") {
    auto selector = f.fetch("chr1:750-780");
    const auto pixels = selector.read_all<T>();
    REQUIRE(pixels.size() == 6);

    CHECK(pixels[0].count == 5476);
    CHECK(pixels[1].count == 5477);
    CHECK(pixels[2].count == 5478);
    CHECK(pixels[3].count == 5511);
    CHECK(pixels[4].count == 5512);
    CHECK(pixels[5].count == 5545);
  }

  SECTION("query does not line up with bins") {
    auto selector = f.fetch("chr1:901-927");
    const auto pixels = selector.read_all<T>();
    REQUIRE(pixels.size() == 6);

    CHECK(pixels[0].count == 5896);
    CHECK(pixels[1].count == 5897);
    CHECK(pixels[2].count == 5898);
    CHECK(pixels[3].count == 5916);
    CHECK(pixels[4].count == 5917);
    CHECK(pixels[5].count == 5935);
  }

  SECTION("large query") {
    auto selector = f.fetch("chr1:75-975");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 4186);

    const auto sum = std::accumulate(
        selector.begin<T>(), selector.end<T>(), T(0),
        [&](T accumulator, const ThinPixel<T>& pixel) { return accumulator + pixel.count; });

    CHECK(sum == 13'405'665);
  }

  SECTION("query spans 1 bin") {
    auto selector = f.fetch("chr1:0-9");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 1);
    CHECK(selector.begin<T>()->count == 1);

    selector = f.fetch("chr1:5-7");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 1);
    CHECK(selector.begin<T>()->count == 1);

    selector = f.fetch("chr1:991-1000");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 1);
    CHECK(selector.begin<T>()->count == 6040);

    selector = f.fetch("chr2:50-60");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 1);
    CHECK(selector.begin<T>()->count == 6091);
  }

  SECTION("query spans 1bp") {
    auto selector = f.fetch("chr1:0-1");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 1);
    CHECK(selector.begin<T>()->count == 1);

    selector = f.fetch("chr2:0-1");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 1);
    CHECK(selector.begin<T>()->count == 6051);

    selector = f.fetch("chr1:12-13");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 1);
    CHECK(selector.begin<T>()->count == 111);

    selector = f.fetch("chr1:999-1000");
    REQUIRE(std::distance(selector.begin<T>(), selector.end<T>()) == 1);
    CHECK(selector.begin<T>()->count == 6040);
  }

  SECTION("query spans entire chromosome") {
    auto selector = f.fetch("chr1");

    CHECK(std::distance(selector.begin<T>(), selector.end<T>()) == 5050);
    auto sum = std::accumulate(
        selector.begin<T>(), selector.end<T>(), T(0),
        [&](T accumulator, const ThinPixel<T>& pixel) { return accumulator + pixel.count; });
    CHECK(sum == 14'420'275);

    selector = f.fetch("chr2");

    CHECK(std::distance(selector.begin<T>(), selector.end<T>()) == 55);
    sum = std::accumulate(
        selector.begin<T>(), selector.end<T>(), T(0),
        [&](T accumulator, const ThinPixel<T>& pixel) { return accumulator + pixel.count; });
    CHECK(sum == 334'290);
  }

  SECTION("equality operator") {
    CHECK(f.fetch("chr1:0-1000") == f.fetch("chr1:0-1000"));
    CHECK(f.fetch("chr1:10-1000") != f.fetch("chr1:0-1000"));
  }

  SECTION("overloads return identical results") {
    CHECK(f.fetch("chr1:0-1000") == f.fetch("chr1", 0, 1000));
    CHECK(f.fetch("chr1\t0\t1000", nullptr, File::QUERY_TYPE::BED) == f.fetch("chr1", 0, 1000));
    CHECK(f.fetch("chr1:0-1000", "chr1:0-1000") == f.fetch("chr1", 0, 1000));
    CHECK(f.fetch("chr1\t0\t1000", "chr2\t0\t99", nullptr, File::QUERY_TYPE::BED) ==
          f.fetch("chr1", 0, 1000, "chr2", 0, 99));
    CHECK(f.fetch(0, 100) == f.fetch("chr1", 0, 1000));
    CHECK(f.fetch(0, 100, 20, 30) == f.fetch("chr1", 0, 1000, "chr1", 200, 300));
  }

  SECTION("invalid queries") {
    CHECK_THROWS_WITH(f.fetch(""), Catch::Matchers::Equals("query is empty"));
    CHECK_THROWS_WITH(f.fetch("chr3"), Catch::Matchers::ContainsSubstring("invalid chromosome"));

    CHECK_THROWS_WITH(f.fetch(":0-1"), Catch::Matchers::ContainsSubstring("invalid chromosome"));
    CHECK_THROWS_WITH(f.fetch("-:0-1"), Catch::Matchers::ContainsSubstring("invalid chromosome"));
    CHECK_THROWS_WITH(f.fetch("::0-1"), Catch::Matchers::ContainsSubstring("invalid chromosome"));

    CHECK_THROWS_WITH(f.fetch("chr1:"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch("chr1-"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch("chr1:-"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch("chr1-0-1"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch("chr1:0:1"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch("chr1:01"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch("chr1:-01"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch("chr1:-01"), Catch::Matchers::ContainsSubstring("malformed"));

    CHECK_THROWS_WITH(f.fetch("chr1:-1"),
                      Catch::Matchers::ContainsSubstring("missing start position"));
    CHECK_THROWS_WITH(f.fetch("chr1:0-"),
                      Catch::Matchers::ContainsSubstring("missing end position"));

    CHECK_THROWS_WITH(f.fetch("chr1:4294967296-0"),
                      Catch::Matchers::ContainsSubstring("invalid start position"));
    CHECK_THROWS_WITH(f.fetch("chr1:0-4294967296"),
                      Catch::Matchers::ContainsSubstring("invalid end position"));

    CHECK_THROWS_WITH(f.fetch("chr1:0-0"),
                      Catch::Matchers::ContainsSubstring(
                          "end position should be greater than the start position"));
    CHECK_THROWS_WITH(f.fetch("chr1:10-5"),
                      Catch::Matchers::ContainsSubstring(
                          "end position should be greater than the start position"));
  }
}

}  // namespace hictk::cooler::test::pixel_selector
