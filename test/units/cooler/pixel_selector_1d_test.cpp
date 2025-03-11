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
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::pixel_selector {

static const auto& datadir = hictk::test::datadir;
static const auto& testdir = hictk::test::testdir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
template <typename N>
static std::pair<std::ptrdiff_t, N> generate_test_data(const std::filesystem::path& path,
                                                       const Reference& chroms,
                                                       std::uint32_t bin_size) {
  auto f = File::create<N>(path.string(), chroms, bin_size, true);

  const auto num_bins = f.bins().size();

  std::vector<ThinPixel<N>> pixels;

  N n = 1;
  N sum = 0;
  for (std::uint64_t bin1_id = 0; bin1_id < num_bins; ++bin1_id) {
    for (std::uint64_t bin2_id = bin1_id; bin2_id < num_bins; ++bin2_id) {
      sum += n;
      pixels.emplace_back(ThinPixel<N>{bin1_id, bin2_id, n++});
    }
  }
  f.append_pixels(pixels.begin(), pixels.end());
  return std::make_pair(static_cast<std::ptrdiff_t>(pixels.size()), sum);
}

TEST_CASE("Cooler (fixed bin size): pixel selector 1D queries", "[pixel_selector][short]") {
  const auto path1 = testdir() / "pixel_selector_devel.cool";

  const Reference chroms{Chromosome{0, "chr1", 1000}, Chromosome{1, "chr2", 100}};
  constexpr std::uint32_t bin_size = 10;
  using T = std::uint32_t;

  const auto [expected_nnz, expected_sum] = generate_test_data<T>(path1, chroms, bin_size);

  const File f(path1.string());
  REQUIRE(std::distance(f.begin<T>(), f.end<T>()) == expected_nnz);
  REQUIRE(std::holds_alternative<BinTableFixed>(f.bins().get()));

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

  SECTION("query spans entire genome") {
    auto selector = f.fetch();
    CHECK(std::distance(selector.begin<T>(), selector.end<T>()) == expected_nnz);
    const auto sum = std::accumulate(
        selector.begin<T>(), selector.end<T>(), T(0),
        [&](T accumulator, const ThinPixel<T>& pixel) { return accumulator + pixel.count; });
    CHECK(sum == expected_sum);
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

TEST_CASE("Cooler (variable bin size): pixel selector 1D queries", "[pixel_selector][short]") {
  const auto path1 = datadir / "cooler" / "cooler_variable_bins_test_file.cool";
  using T = std::uint32_t;

  const File f(path1.string());
  REQUIRE(std::holds_alternative<BinTableVariable<>>(f.bins().get()));

  SECTION("query overlaps chrom start") {
    auto selector = f.fetch("chr1:0-20");
    const auto pixels = selector.read_all<T>();
    REQUIRE(pixels.size() == 1);
    const auto& p = pixels.front();

    CHECK(p.coords.bin1.id() == 0);
    CHECK(p.coords.bin2.id() == 2);
    CHECK(p.count == 7);
  }

  SECTION("query overlaps chrom end") {
    auto selector = f.fetch("chr1:20-32");
    const auto pixels = selector.read_all<T>();
    REQUIRE(pixels.size() == 1);
    const auto& p = pixels.front();

    CHECK(p.coords.bin1.id() == 2);
    CHECK(p.coords.bin2.id() == 3);
    CHECK(p.count == 1);
  }

  SECTION("query does not overlap chrom boundaries") {
    auto selector = f.fetch("chr1:15-23");
    const auto pixels = selector.read_all<T>();
    REQUIRE(pixels.empty());
  }

  SECTION("query does not line up with bins") {
    auto selector = f.fetch("chr1:17-27");
    const auto pixels = selector.read_all<T>();
    REQUIRE(pixels.size() == 1);

    const auto& p = pixels.front();

    CHECK(p.coords.bin1.id() == 2);
    CHECK(p.coords.bin2.id() == 3);
    CHECK(p.count == 1);
  }

  SECTION("query spans 1 bin") {
    auto selector = f.fetch("chr1:0-8");
    CHECK(selector.empty());
  }

  SECTION("query spans 1bp") {
    auto selector = f.fetch("chr1:0-1");
    CHECK(selector.empty());
  }

  SECTION("query spans entire chromosome") {
    auto selector = f.fetch("chr1");
    auto pixels = selector.read_all<T>();

    REQUIRE(pixels.size() == 4);
    CHECK(pixels[0].count == 7);
    CHECK(pixels[1].count == 1);
    CHECK(pixels[2].count == 7);
    CHECK(pixels[3].count == 1);

    selector = f.fetch("chr2");
    pixels = selector.read_all<T>();

    REQUIRE(pixels.size() == 3);
    CHECK(pixels[0].count == 5);
    CHECK(pixels[1].count == 5);
    CHECK(pixels[2].count == 6);
  }

  SECTION("query spans entire genome") {
    constexpr std::ptrdiff_t expected_nnz = 19;
    constexpr T expected_sum = 96;
    auto selector = f.fetch();
    CHECK(std::distance(selector.begin<T>(), selector.end<T>()) == expected_nnz);
    const auto sum = std::accumulate(
        selector.begin<T>(), selector.end<T>(), T(0),
        [&](T accumulator, const ThinPixel<T>& pixel) { return accumulator + pixel.count; });
    CHECK(sum == expected_sum);
  }
}

TEST_CASE("Cooler (storage-mode=square): pixel selector 1D queries", "[pixel_selector][short]") {
  const auto path =
      datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/1000";
  using T = std::uint32_t;

  const File f(path.string());

  SECTION("valid queries") {
    const auto sel = f.fetch();
    const auto sum = std::accumulate(
        sel.template begin<T>(), sel.template end<T>(), std::uint64_t{0},
        [&](std::uint64_t accumulator, const auto& p) { return accumulator + p.count; });
    const auto nnz =
        static_cast<std::size_t>(std::distance(sel.template begin<T>(), sel.template end<T>()));
    CHECK(sum == 594'006'205);
    CHECK(nnz == 4'241'909);
  }

  SECTION("invalid queries") {
    CHECK_THROWS(f.fetch("chr1"));
    CHECK_THROWS(f.fetch("chr1", "chr2"));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::pixel_selector
