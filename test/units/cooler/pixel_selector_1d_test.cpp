// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>

#include "hictk/cooler.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::pixel_selector {

template <typename N>
static std::ptrdiff_t generate_test_data(const std::filesystem::path& path, const Reference& chroms,
                                         std::uint32_t bin_size) {
  auto f = File::create_new_cooler<N>(path.string(), chroms, bin_size, true);

  const auto num_bins = f.bins().size();

  std::vector<Pixel<N>> pixels;

  N n = 1;
  for (std::uint64_t bin1_id = 0; bin1_id < num_bins; ++bin1_id) {
    for (std::uint64_t bin2_id = bin1_id; bin2_id < num_bins; ++bin2_id) {
      pixels.emplace_back(Pixel<N>{f.bins(), bin1_id, bin2_id, n++});
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

  auto f = File::open_read_only(path1.string());
  REQUIRE(std::distance(f.begin<T>(), f.end<T>()) == expected_nnz);

  SECTION("query overlaps chrom start") {
    auto selector = f.fetch<T>("chr1:0-20");
    const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(pixels.size() == 3);

    CHECK(pixels[0].count == 1);
    CHECK(pixels[1].count == 2);
    CHECK(pixels[2].count == 111);
  }

  SECTION("query overlaps chrom end") {
    auto selector = f.fetch<T>("chr1:980-1000");
    const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(pixels.size() == 3);

    CHECK(pixels[0].count == 6028);
    CHECK(pixels[1].count == 6029);
    CHECK(pixels[2].count == 6040);
  }

  SECTION("query does not overlap chrom boundaries") {
    auto selector = f.fetch<T>("chr1:750-780");
    const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(pixels.size() == 6);

    CHECK(pixels[0].count == 5476);
    CHECK(pixels[1].count == 5477);
    CHECK(pixels[2].count == 5478);
    CHECK(pixels[3].count == 5511);
    CHECK(pixels[4].count == 5512);
    CHECK(pixels[5].count == 5545);
  }

  SECTION("query does not line up with bins") {
    auto selector = f.fetch<T>("chr1:901-927");
    const std::vector<Pixel<T>> pixels(selector.begin(), selector.end());
    REQUIRE(pixels.size() == 6);

    CHECK(pixels[0].count == 5896);
    CHECK(pixels[1].count == 5897);
    CHECK(pixels[2].count == 5898);
    CHECK(pixels[3].count == 5916);
    CHECK(pixels[4].count == 5917);
    CHECK(pixels[5].count == 5935);
  }

  SECTION("large query") {
    auto selector = f.fetch<T>("chr1:75-975");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 4186);

    const auto sum = std::accumulate(
        selector.begin(), selector.end(), T(0),
        [&](T accumulator, const Pixel<T>& pixel) { return accumulator + pixel.count; });

    CHECK(sum == 13'405'665);
  }

  SECTION("query spans 1 bin") {
    auto selector = f.fetch<T>("chr1:0-9");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 1);

    selector = f.fetch<T>("chr1:5-7");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 1);

    selector = f.fetch<T>("chr1:991-1000");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 6040);

    selector = f.fetch<T>("chr2:50-60");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 6091);
  }

  SECTION("query spans 1bp") {
    auto selector = f.fetch<T>("chr1:0-1");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 1);

    selector = f.fetch<T>("chr2:0-1");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 6051);

    selector = f.fetch<T>("chr1:12-13");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 111);

    selector = f.fetch<T>("chr1:999-1000");
    REQUIRE(std::distance(selector.begin(), selector.end()) == 1);
    CHECK((*selector.begin()).count == 6040);
  }

  SECTION("query spans entire chromosome") {
    auto selector = f.fetch<T>("chr1");

    CHECK(std::distance(selector.begin(), selector.end()) == 5050);
    auto sum = std::accumulate(
        selector.begin(), selector.end(), T(0),
        [&](T accumulator, const Pixel<T>& pixel) { return accumulator + pixel.count; });
    CHECK(sum == 14'420'275);

    selector = f.fetch<T>("chr2");

    CHECK(std::distance(selector.begin(), selector.end()) == 55);
    sum = std::accumulate(
        selector.begin(), selector.end(), T(0),
        [&](T accumulator, const Pixel<T>& pixel) { return accumulator + pixel.count; });
    CHECK(sum == 334'290);
  }

  SECTION("equality operator") {
    CHECK(f.fetch<T>("chr1:0-1000") == f.fetch<T>("chr1:0-1000"));
    CHECK(f.fetch<T>("chr1:10-1000") != f.fetch<T>("chr1:0-1000"));
  }

  SECTION("overloads return identical results") {
    CHECK(f.fetch<T>("chr1:0-1000") == f.fetch<T>("chr1", 0, 1000));
    CHECK(f.fetch<T>("chr1\t0\t1000", File::QUERY_TYPE::BED) == f.fetch<T>("chr1", 0, 1000));
    CHECK(f.fetch<T>("chr1:0-1000", "chr1:0-1000") == f.fetch<T>("chr1", 0, 1000));
    CHECK(f.fetch<T>("chr1\t0\t1000", "chr2\t0\t99", File::QUERY_TYPE::BED) ==
          f.fetch<T>("chr1", 0, 1000, "chr2", 0, 99));
  }

  SECTION("invalid queries") {
    CHECK_THROWS_WITH(f.fetch<T>(""), Catch::Matchers::Equals("query is empty"));
    CHECK_THROWS_WITH(f.fetch<T>("chr3"), Catch::Matchers::ContainsSubstring("invalid chromosome"));

    CHECK_THROWS_WITH(f.fetch<T>(":0-1"), Catch::Matchers::ContainsSubstring("invalid chromosome"));
    CHECK_THROWS_WITH(f.fetch<T>("-:0-1"),
                      Catch::Matchers::ContainsSubstring("invalid chromosome"));
    CHECK_THROWS_WITH(f.fetch<T>("::0-1"),
                      Catch::Matchers::ContainsSubstring("invalid chromosome"));

    CHECK_THROWS_WITH(f.fetch<T>("chr1:"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1-"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:-"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1-0-1"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:0:1"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:01"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:-01"), Catch::Matchers::ContainsSubstring("malformed"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:-01"), Catch::Matchers::ContainsSubstring("malformed"));

    CHECK_THROWS_WITH(f.fetch<T>("chr1:-1"),
                      Catch::Matchers::ContainsSubstring("missing start position"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:0-"),
                      Catch::Matchers::ContainsSubstring("missing end position"));

    CHECK_THROWS_WITH(f.fetch<T>("chr1:4294967296-0"),
                      Catch::Matchers::ContainsSubstring("invalid start position"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:0-4294967296"),
                      Catch::Matchers::ContainsSubstring("invalid end position"));

    CHECK_THROWS_WITH(f.fetch<T>("chr1:0-0"),
                      Catch::Matchers::ContainsSubstring(
                          "end position should be greater than the start position"));
    CHECK_THROWS_WITH(f.fetch<T>("chr1:10-5"),
                      Catch::Matchers::ContainsSubstring(
                          "end position should be greater than the start position"));
  }
}
}  // namespace hictk::cooler::test::pixel_selector
