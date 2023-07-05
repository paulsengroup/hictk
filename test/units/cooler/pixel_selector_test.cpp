// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/pixel_selector.hpp"

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>

#include "hictk/cooler.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk::test {
inline const internal::TmpDir testdir{true};                     // NOLINT(cert-err58-cpp)
inline const std::filesystem::path datadir{"test/data/cooler"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::cooler::test::pixel_selector {
const auto& testdir = hictk::test::testdir;
const auto& datadir = hictk::test::datadir;

template <typename N>
static std::ptrdiff_t generate_test_data(const std::filesystem::path& path, const Reference& chroms,
                                         std::uint32_t bin_size) {
  auto f = File::create_new_cooler<N>(path.string(), chroms, bin_size, true);

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

  auto f = File::open_read_only(path1.string());
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
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: pixel selector 2D queries", "[pixel_selector][short]") {
  using T = std::uint32_t;
  const auto path = datadir / "cooler_test_file.cool";
  auto f = File::open_read_only(path.string());

  SECTION("cis") {
    SECTION("overloads return identical results") {
      CHECK(f.fetch("1:5000000-5500000", "1:5000000-6500000") ==
            f.fetch("1", 5000000, 5500000, "1", 5000000, 6500000));
    }

    SECTION("valid") {
      auto selector = f.fetch("1:5000000-5500000", "1:5000000-6500000");
      const auto pixels = selector.read_all<T>();
      REQUIRE(pixels.size() == 8);

      CHECK(pixels[0].count == 20);
      CHECK(pixels[1].count == 1);
      CHECK(pixels[2].count == 18);
      CHECK(pixels[3].count == 8);
      CHECK(pixels[4].count == 1);
      CHECK(pixels[5].count == 9);
      CHECK(pixels[6].count == 6);
      CHECK(pixels[7].count == 2);
    }

    SECTION("empty") {
      auto selector = f.fetch("1:0-100000");
      CHECK(selector.begin<T>() == selector.end<T>());
    }
  }

  SECTION("trans") {
    SECTION("overloads return identical results") {
      CHECK(f.fetch("1:48000000-50000000", "4:30000000-35000000") ==
            f.fetch("1", 48000000, 50000000, "4", 30000000, 35000000));
    }
    SECTION("valid") {
      auto selector = f.fetch("1:48000000-50000000", "4:30000000-35000000");
      const auto pixels = selector.read_all<T>();
      REQUIRE(pixels.size() == 6);

      CHECK(pixels[0].count == 1);
      CHECK(pixels[1].count == 3);
      CHECK(pixels[2].count == 1);
      CHECK(pixels[3].count == 3);
      CHECK(pixels[4].count == 7);
      CHECK(pixels[5].count == 1);
    }

    SECTION("empty") {
      auto selector = f.fetch("1:0-50000", "2:0-50000");
      CHECK(selector.begin<T>() == selector.end<T>());
    }
  }
}
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: pixel selector w/ balancing", "[pixel_selector][short]") {
  auto path = datadir / "ENCFF993FGR.2500000.cool";
  auto clr = File::open_read_only(path.string());

  SECTION("read weights") {
    SECTION("valid") {
      CHECK(clr.read_weights("weight")->type() == balancing::Weights::Type::MULTIPLICATIVE);
      for (const auto* name : {"GW_SCALE", "INTER_SCALE", "SCALE", "VC", "VC_SQRT"}) {
        CHECK(clr.read_weights(name)->type() == balancing::Weights::Type::DIVISIVE);
      }
    }

    SECTION("invalid") {
      CHECK_THROWS(clr.read_weights(""));
      CHECK_THROWS(clr.read_weights("AAA"));
    }

    SECTION("purging") {
      CHECK(clr.purge_weights() == false);
      CHECK(clr.purge_weights("weight") == false);

      const auto w = clr.read_weights("weight");
      CHECK(w.use_count() == 2);
      CHECK(clr.purge_weights("weight") == true);
      CHECK(w.use_count() == 1);

      clr.read_weights("weight");
      CHECK(clr.purge_weights() == true);
    }
  }

  SECTION("1D query") {
    const auto selector = clr.fetch("chr1", 5'000'000, 10'000'000, clr.read_weights("weight"));
    constexpr std::array<double, 3> expected{3.345797, 0.328794, 4.456354};
    const auto pixels = selector.read_all<double>();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK_THAT(pixels[i].count, Catch::Matchers::WithinAbs(expected[i], 1.0e-6));
    }
  }

  SECTION("2D query") {
    const auto selector = clr.fetch("chr1", 5'000'000, 10'000'000, "chr2", 5'000'000, 10'000'000,
                                    clr.read_weights("weight"));
    constexpr std::array<double, 4> expected{0.001782, 0.002756, 0.002047, 0.004749};
    const auto pixels = selector.read_all<double>();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK_THAT(pixels[i].count, Catch::Matchers::WithinAbs(expected[i], 1.0e-6));
    }
  }

  SECTION("invalid iterator type") {
    const auto selector = clr.fetch("chr1", 5'000'000, 10'000'000, "chr2", 5'000'000, 10'000'000,
                                    clr.read_weights("weight"));
    CHECK_THROWS(selector.read_all<std::int32_t>());
  }
}

}  // namespace hictk::cooler::test::pixel_selector
