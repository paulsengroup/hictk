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

const auto pathV8 =
    (hictk::test::datadir / "4DNFIZ1ZVXC8.hic8").string();  // NOLINT(cert-err58-cpp)
const auto pathV9 =
    (hictk::test::datadir / "4DNFIZ1ZVXC8.hic9").string();              // NOLINT(cert-err58-cpp)
const auto path_binary = (hictk::test::datadir / "data.zip").string();  // NOLINT(cert-err58-cpp)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static std::vector<hictk::Pixel<float>> head(const std::vector<hictk::Pixel<float>>& buffer,
                                             std::size_t n = 5) {
  REQUIRE(buffer.size() >= n);

  std::vector<hictk::Pixel<float>> slice(n);
  std::copy_n(buffer.begin(), n, slice.begin());
  return slice;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static std::vector<hictk::Pixel<float>> tail(const std::vector<hictk::Pixel<float>>& buffer,
                                             std::size_t n = 5) {
  REQUIRE(buffer.size() >= n);

  std::vector<hictk::Pixel<float>> slice(n);
  std::copy_n(buffer.end() - std::int32_t(n), n, slice.begin());
  return slice;
}

template <typename N>
static N sumCounts(const std::vector<hictk::Pixel<float>>& buffer) {
  return std::accumulate(buffer.begin(), buffer.end(), N(0),
                         [](N accumulator, const hictk::Pixel<float>& p) {
                           return accumulator + static_cast<N>(p.count);
                         });
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void compareContactRecord(const hictk::Pixel<float>& r1, const SerializedPixel& r2) {
  CHECK(r1.coords.bin1.start() == r2.bin1_id);
  CHECK(r1.coords.bin2.start() == r2.bin2_id);
  CHECK_THAT(r1.count, Catch::Matchers::WithinRel(r2.count));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (observed NONE BP 10000)", "[hic][short]") {
  SECTION("intra-chromosomal") {
    constexpr std::size_t expected_size = 1433133;
    constexpr std::int32_t expected_sum = 19968156;

    constexpr std::size_t N = 5;
    constexpr std::array<float, N> head_expected{1745, 2844, 409, 195, 195};
    constexpr std::array<float, N> tail_expected{119, 34, 281, 53, 193};

    constexpr auto expected_value =
        std::make_pair(std::size_t(1229799), SerializedPixel{15770000, 15770000, 1234.0F});

    SECTION("v8") {
      auto sel = HiCFile(pathV8, 10'000, MatrixType::observed, MatrixUnit::BP).fetch("chr2L");
      const auto buffer = sel.read_all<float>();
      REQUIRE(buffer.size() == expected_size);

      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
        CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }

    SECTION("v9") {
      auto sel = HiCFile(pathV9, 10'000, MatrixType::observed, MatrixUnit::BP).fetch("chr2L");
      const auto buffer = sel.read_all<float>();
      REQUIRE(buffer.size() == expected_size);

      CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

      const auto h = head(buffer, N);
      const auto t = tail(buffer, N);

      for (std::size_t i = 0; i < N; ++i) {
        CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
        CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
      }

      compareContactRecord(buffer[expected_value.first], expected_value.second);
      CHECK(std::is_sorted(buffer.begin(), buffer.end()));
    }
  }
}
