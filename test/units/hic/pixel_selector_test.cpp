// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <filesystem>
#include <string>

#include "hictk/hic.hpp"

using namespace hictk::hic;

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

template <typename N>
using Pixel = hictk::Pixel<N>;

// NOLINTNEXTLINE(cert-err58-cpp)
const auto pathV8 = (hictk::test::datadir / "4DNFIZ1ZVXC8.hic8").string();

// NOLINTNEXTLINE(cert-err58-cpp)
const auto pathV9 = (hictk::test::datadir / "4DNFIZ1ZVXC8.hic9").string();

// NOLINTNEXTLINE(cert-err58-cpp)
const auto path_binary = (hictk::test::datadir / "data.zip").string();

template <typename N>
static std::vector<hictk::Pixel<N>> head(const std::vector<hictk::Pixel<N>>& buffer,
                                         std::size_t n = 5) {
  REQUIRE(buffer.size() >= n);

  std::vector<hictk::Pixel<N>> slice(n);
  std::copy_n(buffer.begin(), n, slice.begin());
  return slice;
}

template <typename N>  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
static std::vector<hictk::Pixel<N>> tail(const std::vector<hictk::Pixel<N>>& buffer,
                                         std::size_t n = 5) {
  REQUIRE(buffer.size() >= n);

  std::vector<hictk::Pixel<N>> slice(n);
  std::copy_n(buffer.end() - std::int32_t(n), n, slice.begin());
  return slice;
}

template <typename N>
static N sumCounts(const std::vector<hictk::Pixel<N>>& buffer) {
  return std::accumulate(buffer.begin(), buffer.end(), N(0),
                         [](N accumulator, const hictk::Pixel<N>& p) {
                           return accumulator + static_cast<N>(p.count);
                         });
}

template <typename N>  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void checkContactRecordsAreWithinBound(std::uint32_t start1, std::uint32_t end1,
                                              std::uint32_t start2, std::uint32_t end2,
                                              const std::vector<Pixel<N>>& buffer) {
  assert(start1 < end1);
  assert(start2 < end2);

  for (const auto& r : buffer) {
    CHECK(r.coords.bin1.start() >= std::min(start1, start2));
    CHECK(r.coords.bin1.end() < std::max(end1, end2));
    CHECK(r.coords.bin2.start() >= std::min(start1, start2));
    CHECK(r.coords.bin2.end() < std::max(end1, end2));
  }
}

template <typename N>  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void compareContactRecord(const hictk::Pixel<N>& r1, const SerializedPixel& r2) {
  CHECK(r1.coords.bin1.start() == r2.bin1_id);
  CHECK(r1.coords.bin2.start() == r2.bin2_id);
  if constexpr (std::is_floating_point_v<N>) {
    CHECK_THAT(r1.count, Catch::Matchers::WithinRel(r2.count));
  } else {
    CHECK(r1.count == static_cast<N>(r2.count));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: pixel selector accessors", "[hic][short]") {
  const auto sel = HiCFile(pathV8, 2'500'000, MatrixType::observed, MatrixUnit::BP)
                       .fetch("chr2L", NormalizationMethod::NONE);

  CHECK(sel.chrom1().name() == "chr2L");
  CHECK(sel.chrom2().name() == "chr2L");
  CHECK(sel.matrix_type() == MatrixType::observed);
  CHECK(sel.normalization() == NormalizationMethod::NONE);
  CHECK(sel.unit() == MatrixUnit::BP);
  CHECK(sel.resolution() == 2500000);

  REQUIRE(sel.chrom1().size() == 23513712);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: pixel selector fetch (observed NONE BP 10000)", "[hic][long]") {
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr std::int32_t expected_sum = 19968156;

    constexpr std::size_t N = 5;
    constexpr std::array<std::int32_t, N> head_expected{1745, 2844, 409, 195, 195};
    constexpr std::array<std::int32_t, N> tail_expected{119, 34, 281, 53, 193};

    constexpr auto expected_value =
        std::make_pair(std::size_t(1229799), SerializedPixel{15770000, 15770000, 1234.0F});

    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP).fetch("chr2L");
      const auto buffer = sel.read_all<std::int32_t>();
      REQUIRE(buffer.size() == expected_size);

      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK(head_expected[i] == h[i].count);
        CHECK(tail_expected[i] == t[i].count);
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP).fetch("chr2L");
      const auto buffer = sel.read_all<std::int32_t>();
      REQUIRE(buffer.size() == expected_size);

      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK(head_expected[i] == h[i].count);
        CHECK(tail_expected[i] == t[i].count);
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }
  }

  SECTION("inter-chromosomal") {
    constexpr std::size_t expected_size = 56743;
    constexpr std::int32_t expected_sum = 70567;

    constexpr std::size_t N = 5;
    constexpr std::array<std::int32_t, N> head_expected{1, 1, 1, 1, 1};
    constexpr std::array<std::int32_t, N> tail_expected{1, 1, 1, 1, 1};

    constexpr auto expected_value =
        std::make_pair(std::size_t(3541), SerializedPixel{770000, 1300000, 13.0F});

    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L", "chr4", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<std::int32_t>();
      REQUIRE(buffer.size() == expected_size);

      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK(head_expected[i] == h[i].count);
        CHECK(tail_expected[i] == t[i].count);
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L", "chr4", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<std::int32_t>();
      REQUIRE(buffer.size() == expected_size);

      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK(head_expected[i] == h[i].count);
        CHECK(tail_expected[i] == t[i].count);
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }
  }

  SECTION("cover type 2 interactions") {
    auto sel = HiCFile(pathV8, 2'500'000, MatrixType::observed, MatrixUnit::BP)
                   .fetch("chr2L", "chr2R", NormalizationMethod::NONE);
    const auto buffer = sel.read_all<std::int32_t>();
    REQUIRE(buffer.size() == 110);
    CHECK(sumCounts<std::int32_t>(buffer) == 1483112);

    compareContactRecord(buffer[38], SerializedPixel{7500000, 12500000, 16512});
    CHECK(std::is_sorted(buffer.begin(), buffer.end()));
  }

  SECTION("sub-chromosomal queries") {
    const std::uint32_t resolution = 10'000;
    SECTION("single pixel") {
      auto sel = HiCFile(pathV9, resolution, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L:100,000-100,001", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<std::int32_t>();
      REQUIRE(buffer.size() == 1);
      compareContactRecord(buffer.front(), SerializedPixel{100000, 100000, 13895.0F});
    }

    SECTION("upper-triangle") {
      auto sel = HiCFile(pathV9, resolution, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L:123,456-200,000", "chr2L:0-200,000", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<std::int32_t>();
      REQUIRE(buffer.size() == 36);
      CHECK(sumCounts<std::int32_t>(buffer) == 99946);
      compareContactRecord(buffer[33], SerializedPixel{180000, 180000, 3888});

      checkContactRecordsAreWithinBound(123456, 200000 + resolution, 0, 200000 + resolution,
                                        buffer);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }

    SECTION("lower-triangle") {
      auto sel = HiCFile(pathV9, resolution, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L:0-200,000", "chr2L:123,456-200,000", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<std::int32_t>();
      REQUIRE(buffer.size() == 132);
      CHECK(sumCounts<std::int32_t>(buffer) == 124561);
      compareContactRecord(buffer[33], SerializedPixel{40000, 130000, 148});
      checkContactRecordsAreWithinBound(0, 200000 + resolution, 123456, 200000 + resolution,
                                        buffer);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }

    SECTION("inter-chromosomal") {
      auto sel = HiCFile(pathV9, resolution, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L:123,456-200,000", "chr4:0-200,000", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<std::int32_t>();
      REQUIRE(buffer.size() == 57);
      CHECK(sumCounts<std::int32_t>(buffer) == 74);
      checkContactRecordsAreWithinBound(123456, 200000 + resolution, 0, 200000 + resolution,
                                        buffer);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }
  }
  SECTION("invalid") {
    SECTION("invalid chromosome") {
      const HiCFile hic(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP);
      CHECK_THROWS(hic.fetch("chr123", NormalizationMethod::NONE));
    }
    SECTION("invalid unit") {
      const HiCFile hic(pathV9, 10'000, MatrixType::observed, MatrixUnit::FRAG);
      CHECK_THROWS(hic.fetch("chr2L", NormalizationMethod::NONE));
    }
    SECTION("expected + norm") {
      const HiCFile hic(pathV9, 10'000, MatrixType::expected, MatrixUnit::BP);
      CHECK_THROWS(hic.fetch("chr2L", NormalizationMethod::VC));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: pixel selector fetch (observed VC BP 10000)", "[hic][long]") {
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr double expected_sum = 20391277.41514;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L", NormalizationMethod::VC);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L", NormalizationMethod::VC);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
  SECTION("inter-chromosomal") {
    constexpr std::size_t expected_size = 56743;
    constexpr double expected_sum = 96690.056244753;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L", "chr4", NormalizationMethod::VC);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP)
                     .fetch("chr2L", "chr4", NormalizationMethod::VC);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: pixel selector fetch (expected NONE BP 10000)", "[hic][long]") {
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr double expected_sum = 18314748.068024;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::expected, MatrixUnit::BP)
                     .fetch("chr2L", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::expected, MatrixUnit::BP)
                     .fetch("chr2L", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
  SECTION("inter-chromosomal") {
    constexpr std::size_t expected_size = 56743;
    constexpr double expected_sum = 12610.80619812;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::expected, MatrixUnit::BP)
                     .fetch("chr2L", "chr4", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::expected, MatrixUnit::BP)
                     .fetch("chr2L", "chr4", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: pixel selector fetch (oe NONE BP 10000)", "[hic][long]") {
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr double expected_sum = 2785506.2274201;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::oe, MatrixUnit::BP)
                     .fetch("chr2L", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::oe, MatrixUnit::BP)
                     .fetch("chr2L", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
  SECTION("inter-chromosomal") {
    constexpr std::size_t expected_size = 56743;
    constexpr double expected_sum = 317520.00459671;
    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::oe, MatrixUnit::BP)
                     .fetch("chr2L", "chr4", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::oe, MatrixUnit::BP)
                     .fetch("chr2L", "chr4", NormalizationMethod::NONE);
      const auto buffer = sel.read_all<double>();
      REQUIRE(buffer.size() == expected_size);
      CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: pixel selector fetch all (observed NONE BP 100000)", "[hic][long]") {
  SECTION("accessors") {
    auto sel = HiCFile(pathV8, 100'000, MatrixType::observed, MatrixUnit::BP).fetch();

    CHECK(sel.resolution() == 100'000);
    CHECK(sel.matrix_type() == MatrixType::observed);
    CHECK(sel.normalization() == NormalizationMethod::NONE);
    CHECK(sel.unit() == MatrixUnit::BP);
    CHECK(sel.bins().size() == 1382);
  }

  SECTION("v8") {
    auto sel = HiCFile(pathV8, 100'000, MatrixType::observed, MatrixUnit::BP).fetch();
    const auto buffer = sel.read_all<double>();
    REQUIRE(buffer.size() == 890384);

    CHECK_THAT(sumCounts(buffer), Catch::Matchers::WithinRel(119208613, 1.0e-6));
    CHECK(std::is_sorted(buffer.begin(), buffer.end()));
  }

  SECTION("v9") {
    auto sel = HiCFile(pathV9, 100'000, MatrixType::observed, MatrixUnit::BP).fetch();
    const auto buffer = sel.read_all<double>();
    REQUIRE(buffer.size() == 890384);

    CHECK_THAT(sumCounts(buffer), Catch::Matchers::WithinRel(119208613, 1.0e-6));
    CHECK(std::is_sorted(buffer.begin(), buffer.end()));
  }
}
