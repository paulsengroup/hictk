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

const auto pathV8 =
    (hictk::test::datadir / "4DNFIZ1ZVXC8.hic8").string();  // NOLINT(cert-err58-cpp)
const auto pathV9 =
    (hictk::test::datadir / "4DNFIZ1ZVXC8.hic9").string();              // NOLINT(cert-err58-cpp)
const auto path_binary = (hictk::test::datadir / "data.zip").string();  // NOLINT(cert-err58-cpp)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
template <typename N>
static std::vector<hictk::Pixel<N>> head(const std::vector<hictk::Pixel<N>>& buffer,
                                         std::size_t n = 5) {
  REQUIRE(buffer.size() >= n);

  std::vector<hictk::Pixel<N>> slice(n);
  std::copy_n(buffer.begin(), n, slice.begin());
  return slice;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
template <typename N>
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
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
template <typename N>
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
template <typename N>
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
TEST_CASE("MatrixSelector accessors", "[hic][short]") {
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
TEST_CASE("MatrixSelector LRU cache", "[hic][short]") {
  HiCFile f(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP);

  auto sel = f.fetch("chr2L", NormalizationMethod::NONE, HiCFile::QUERY_TYPE::UCSC, 100);
  // Fill cache
  const auto expected_sum = sumCounts(sel.read_all<std::int32_t>());

  REQUIRE(f.block_cache_hit_rate() == 0);
  REQUIRE(f.block_cache_size() == 6);

  auto sum = sumCounts<std::int32_t>(sel.read_all<std::int32_t>());
  CHECK(sum == expected_sum);
  CHECK(f.block_cache_hit_rate() == 0.5);
  CHECK(f.block_cache_size() == 6);

  for (auto i = 0; i < 3; ++i) {
    sum = sumCounts<std::int32_t>(sel.read_all<std::int32_t>());
    CHECK(sum == expected_sum);
  }
  CHECK(f.block_cache_hit_rate() == 4.0 / 5.0);
  CHECK(f.block_cache_size() == 6);

  f.clear_block_cache();
  CHECK(f.block_cache_hit_rate() == 0);
  CHECK(f.block_cache_size() == 0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (observed NONE BP 10000)", "[hic][short]") {
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
      auto f = std::fopen("/tmp/test.bg2", "w");
      fmt::print(f, FMT_STRING("{}\n"), fmt::join(buffer, "\n"));
      std::fclose(f);
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
}
