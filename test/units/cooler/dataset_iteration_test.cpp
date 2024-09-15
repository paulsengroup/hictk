// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <iterator>
#include <random>
#include <vector>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::dataset {
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset linear iteration", "[dataset][long]") {
  const auto path = datadir / "cooler_test_file.cool";

  const RootGroup grp{HighFive::File(path.string()).getGroup("/")};
  const Dataset dset(grp, "/pixels/count");

  std::vector<std::uint32_t> pixel_buff;
  dset.read_all(pixel_buff);
  REQUIRE(pixel_buff.size() == 107'041);

  SECTION("forward") {
    auto it = dset.begin<std::uint32_t>(32'000);
    auto last_pixel = dset.end<std::uint32_t>(32'000);
    REQUIRE(std::distance(it, last_pixel) == 107'041);

    for (const auto& expected : pixel_buff) {
      REQUIRE(it != last_pixel);
      CHECK(*it++ == expected);
    }
    CHECK(it == last_pixel);
  }

  SECTION("backward") {
    auto it = dset.end<std::uint32_t>(32'000);
    auto first_pixel = dset.begin<std::uint32_t>(32'000);

    REQUIRE(std::distance(first_pixel, it) == 107'041);

    for (std::size_t i = pixel_buff.size(); i != 0; --i) {
      REQUIRE(it != first_pixel);
      CHECK(*(--it) == pixel_buff[i - 1]);
    }

    CHECK(it == first_pixel);
  }
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset random iteration", "[dataset][medium]") {
  const auto path = testdir() / "dataset_iterator_random.h5";

  const RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
  Dataset dset(grp, "int", std::uint8_t{}, HighFive::DataSpace::UNLIMITED,
               Dataset::default_access_props(), Dataset::init_create_props(1, 32'000));

  std::random_device rd;
  std::mt19937_64 rand_eng{rd()};

  constexpr std::size_t N = 5'000'000;
  constexpr auto M = static_cast<std::ptrdiff_t>(N) / 2;
  std::vector<std::uint8_t> buff(N);
  std::generate(buff.begin(), buff.end(),
                [&]() { return std::uniform_int_distribution<std::uint8_t>{0, 255}(rand_eng); });
  dset.append(buff);
  REQUIRE(dset.size() == N);

  SECTION("operator -/+") {
    const auto first = dset.begin<std::uint8_t>(32'000);
    const auto mid = first + M;
    const auto last = dset.end<std::uint8_t>(32'000);
    for (std::size_t i = 0; i < 10000; ++i) {
      auto js = std::uniform_int_distribution<std::ptrdiff_t>{0, N - 1}(rand_eng);
      auto ju = static_cast<std::size_t>(js);

      CHECK(*(first + js) == buff[ju]);
      CHECK(*(last - js) == buff[N - ju]);

      js = std::uniform_int_distribution<std::ptrdiff_t>{-(M - 1), M - 1}(rand_eng);
      const auto ju1 = static_cast<std::size_t>(M + js);
      const auto ju2 = static_cast<std::size_t>(M - js);
      CHECK(*(mid + js) == buff[ju1]);
      CHECK(*(mid - js) == buff[ju2]);
    }
  }

  SECTION("subsequent calls to operator+=") {
    for (std::size_t i = 0; i < 50; ++i) {
      auto first = dset.begin<std::uint8_t>(32'000);
      auto last = dset.end<std::uint8_t>(32'000);
      std::ptrdiff_t j = 0;

      while (first < last) {
        CHECK(*first == buff[static_cast<std::size_t>(j)]);

        const auto lb = std::max(std::ptrdiff_t{-100}, -j);
        const auto ub = std::min(std::ptrdiff_t{500}, static_cast<std::ptrdiff_t>(buff.size()) - j);
        const auto step = std::uniform_int_distribution<std::ptrdiff_t>{lb, ub}(rand_eng);
        first += step;
        j += step;
      }
    }
  }

  SECTION("subsequent calls to operator-=") {
    for (std::size_t i = 0; i < 50; ++i) {
      auto first = dset.end<std::uint8_t>(32'000) - 1;
      auto last = dset.begin<std::uint8_t>(32'000);
      auto j = static_cast<std::ptrdiff_t>(buff.size() - 1);

      while (first > last) {
        CHECK(*first == buff[static_cast<std::size_t>(j)]);

        const auto lb =
            std::max(std::ptrdiff_t{-100}, -(static_cast<std::ptrdiff_t>(buff.size()) - j));
        const auto ub = std::min(std::ptrdiff_t{500}, j);
        const auto step = std::uniform_int_distribution<std::ptrdiff_t>{lb, ub}(rand_eng);

        first -= step;
        j -= step;
      }
    }
  }
}

}  // namespace hictk::cooler::test::dataset
