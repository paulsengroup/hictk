// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <filesystem>
#include <numeric>
#include <string>

#include "hictk/fmt.hpp"
#include "hictk/hic.hpp"

using namespace hictk;

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

const auto pathV8 = (test::datadir / "4DNFIZ1ZVXC8.hic8").string();  // NOLINT(cert-err58-cpp)
const auto pathV9 = (test::datadir / "4DNFIZ1ZVXC8.hic9").string();  // NOLINT(cert-err58-cpp)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static std::vector<Pixel<float>> head(const std::vector<Pixel<float>>& buffer, std::size_t n = 5) {
  REQUIRE(buffer.size() >= n);

  std::vector<Pixel<float>> slice(n);
  std::copy_n(buffer.begin(), n, slice.begin());
  return slice;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static std::vector<Pixel<float>> tail(const std::vector<Pixel<float>>& buffer, std::size_t n = 5) {
  REQUIRE(buffer.size() >= n);

  std::vector<Pixel<float>> slice(n);
  std::copy_n(buffer.end() - std::int32_t(n), n, slice.begin());
  return slice;
}

template <typename N>
static N sumCounts(const std::vector<Pixel<float>>& buffer) {
  return std::accumulate(
      buffer.begin(), buffer.end(), N(0),
      [](N accumulator, const Pixel<float>& r) { return accumulator + static_cast<N>(r.count); });
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void checkContactRecordsAreWithinBound(std::uint32_t start1, std::uint32_t end1,
                                              std::uint32_t start2, std::uint32_t end2,
                                              const std::vector<Pixel<float>>& buffer) {
  assert(start1 < end1);
  assert(start2 < end2);

  for (const auto& r : buffer) {
    CHECK(r.coords.bin1.start() >= std::min(start1, start2));
    CHECK(r.coords.bin1.end() < std::max(end1, end2));
    CHECK(r.coords.bin2.start() >= std::min(start1, start2));
    CHECK(r.coords.bin2.end() < std::max(end1, end2));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void compareContactRecord(const Pixel<float>& r1, const SerializedPixel& r2) {
  CHECK(r1.coords.bin1.start() == r2.bin1_id);
  CHECK(r1.coords.bin2.start() == r2.bin2_id);
  CHECK_THAT(r1.count, Catch::Matchers::WithinRel(r2.count));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector accessors", "[hic][short]") {
  const auto sel = HiCFile(pathV8, 2'500'000, MatrixType::observed, MatrixUnit::BP)
                       .get_matrix_selector("chr2L", NormalizationMethod::NONE);

  CHECK(sel.chrom1().name() == "chr2L");
  CHECK(sel.chrom2().name() == "chr2L");
  CHECK(sel.matrix_type() == MatrixType::observed);
  CHECK(sel.normalizationMethod() == NormalizationMethod::NONE);
  CHECK(sel.matrixUnit() == MatrixUnit::BP);
  CHECK(sel.resolution() == 2500000);

  REQUIRE(sel.chrom1().size() == 23513712);
  CHECK(sel.numBins1() == 10);
  CHECK(sel.numBins2() == 10);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector LRU cache", "[hic][short]") {
  std::vector<Pixel<float>> buffer;
  HiCFile f(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP);

  auto sel = f.get_matrix_selector("chr2L", NormalizationMethod::NONE);

  CHECK(sel.blockCacheHitRate() == 0.0);
  CHECK(sel.blockCacheSize() == 0);

  // Fill cache
  sel.fetch(buffer);
  CHECK(sel.blockCacheHitRate() == 0.0);

  sel.fetch(buffer);
  CHECK(sel.blockCacheHitRate() == 0.5);
  CHECK(sel.blockCacheSize() == 6);

  for (auto i = 0; i < 5; ++i) {
    sel.fetch(buffer);
  }
  CHECK(sel.blockCacheHitRate() == 6.0 / 7.0);
  CHECK(sel.blockCacheSize() == 6);

  sel.clearBlockCache();
  CHECK(sel.blockCacheHitRate() == 0);
  CHECK(sel.blockCacheSize() == 0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (observed NONE BP 10000)", "[hic][short]") {
  std::vector<Pixel<float>> buffer;
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr std::int32_t expected_sum = 19968156;

    constexpr std::size_t N = 5;
    constexpr std::array<float, N> head_expected{1745, 2844, 409, 195, 195};
    constexpr std::array<float, N> tail_expected{119, 34, 281, 53, 193};

    constexpr auto expected_value =
        std::make_pair(std::size_t(1229799), SerializedPixel{15770000, 15770000, 1234.0F});

    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
        CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
    }
    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
        CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
    }
  }

  SECTION("inter-chromosomal") {
    constexpr std::size_t expected_size = 56743;
    constexpr std::int32_t expected_sum = 70567;

    constexpr std::size_t N = 5;
    constexpr std::array<float, N> head_expected{1, 1, 1, 1, 1};
    constexpr std::array<float, N> tail_expected{1, 1, 1, 1, 1};

    constexpr auto expected_value =
        std::make_pair(std::size_t(3541), SerializedPixel{770000, 1300000, 13.0F});

    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr4", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
        CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr4", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
        CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
    }

    SECTION("cover type 2 interactions") {
      auto sel = HiCFile(pathV8, 2'500'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr2R", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == 110);
      CHECK(sumCounts<std::int32_t>(buffer) == 1483112);

      compareContactRecord(buffer[38], SerializedPixel{7500000, 12500000, 16512});
    }

    SECTION("sub-queries") {
      const std::uint32_t resolution = 10'000;
      SECTION("single pixel") {
        auto sel = HiCFile(pathV9, resolution, MatrixType::observed, MatrixUnit::BP)
                       .get_matrix_selector("chr2L", NormalizationMethod::NONE);
        sel.fetch(100000, 100001, 100000, 100001, buffer);
        REQUIRE(buffer.size() == 1);
        compareContactRecord(buffer.front(), SerializedPixel{100000, 100000, 13895.0F});
      }

      SECTION("upper-triangle") {
        auto sel = HiCFile(pathV9, resolution, MatrixType::observed, MatrixUnit::BP)
                       .get_matrix_selector("chr2L", NormalizationMethod::NONE);
        sel.fetch(123456, 200000, 0, 200000, buffer, true);
        REQUIRE(buffer.size() == 132);
        CHECK(sumCounts<std::int32_t>(buffer) == 124561);
        compareContactRecord(buffer[33], SerializedPixel{40000, 130000, 148});
        checkContactRecordsAreWithinBound(123456, 200000 + resolution, 0, 200000 + resolution,
                                          buffer);
      }

      SECTION("lower-triangle") {
        auto sel = HiCFile(pathV9, resolution, MatrixType::observed, MatrixUnit::BP)
                       .get_matrix_selector("chr2L", NormalizationMethod::NONE);
        sel.fetch(0, 200000, 123456, 200000, buffer, true);
        REQUIRE(buffer.size() == 132);
        CHECK(sumCounts<std::int32_t>(buffer) == 124561);
        compareContactRecord(buffer[33], SerializedPixel{40000, 130000, 148});
        checkContactRecordsAreWithinBound(0, 200000 + resolution, 123456, 200000 + resolution,
                                          buffer);
      }

      SECTION("inter-chromosomal") {
        auto sel = HiCFile(pathV9, resolution, MatrixType::observed, MatrixUnit::BP)
                       .get_matrix_selector("chr2L", "chr4", NormalizationMethod::NONE);
        sel.fetch(123456, 200000, 0, 200000, buffer);
        REQUIRE(buffer.size() == 57);
        CHECK(sumCounts<std::int32_t>(buffer) == 74);
        checkContactRecordsAreWithinBound(123456, 200000 + resolution, 0, 200000 + resolution,
                                          buffer);
      }
    }

    SECTION("invalid") {
      SECTION("invalid chromosome") {
        HiCFile hic(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP);
        CHECK_THROWS(hic.get_matrix_selector("chr123", NormalizationMethod::NONE));
        CHECK_THROWS(hic.get_matrix_selector(999, NormalizationMethod::NONE));
      }
      SECTION("invalid unit") {
        HiCFile hic(pathV9, 10'000, MatrixType::observed, MatrixUnit::FRAG);
        CHECK_THROWS(hic.get_matrix_selector("chr2L", NormalizationMethod::NONE));
      }
      SECTION("expected + norm") {
        HiCFile hic(pathV9, 10'000, MatrixType::expected, MatrixUnit::BP);
        CHECK_THROWS(hic.get_matrix_selector("chr2L", NormalizationMethod::VC));
      }
      SECTION("invalid range") {
        HiCFile hic(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP);
        CHECK_THROWS(
            hic.get_matrix_selector("chr2L", NormalizationMethod::NONE).fetch(1000, 0, buffer));
        CHECK_THROWS(hic.get_matrix_selector("chr2L", NormalizationMethod::NONE)
                         .fetch(0, 1'000'000'000, buffer));
      }
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (observed VC BP 10000)", "[hic][short]") {
  std::vector<Pixel<float>> buffer;
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr double expected_sum = 20391277.41514;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", NormalizationMethod::VC);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", NormalizationMethod::VC);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
  SECTION("inter-chromosomal") {
    constexpr std::size_t expected_size = 56743;
    constexpr double expected_sum = 96690.056244753;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr4", NormalizationMethod::VC);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr4", NormalizationMethod::VC);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (expected NONE BP 10000)", "[hic][short]") {
  std::vector<Pixel<float>> buffer;
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr double expected_sum = 18314748.068024;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::expected, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::expected, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
  SECTION("inter-chromosomal") {
    constexpr std::size_t expected_size = 56743;
    constexpr double expected_sum = 12610.80619812;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::expected, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr4", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::expected, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr4", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (oe NONE BP 10000)", "[hic][short]") {
  std::vector<Pixel<float>> buffer;
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr double expected_sum = 2785506.2274201;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::oe, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::oe, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
  SECTION("inter-chromosomal") {
    constexpr std::size_t expected_size = 56743;
    constexpr double expected_sum = 317520.00459671;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::oe, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr4", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::oe, MatrixUnit::BP)
                     .get_matrix_selector("chr2L", "chr4", NormalizationMethod::NONE);
      sel.fetch(buffer, true);
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
}
