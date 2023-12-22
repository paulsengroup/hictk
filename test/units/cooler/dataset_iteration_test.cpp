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
  Dataset dset(grp, "int", std::uint64_t{});

  std::random_device rd;
  std::mt19937_64 rand_eng{rd()};

  constexpr std::size_t N = 5'000'000;
  std::vector<std::uint64_t> buff(N);
  std::generate(buff.begin(), buff.end(), [&]() { return rand_eng(); });
  dset.append(buff);
  REQUIRE(dset.size() == N);

  SECTION("operator -/+") {
    auto first = dset.begin<std::uint64_t>(32'000);
    auto last = dset.end<std::uint64_t>(32'000);
    for (std::size_t i = 0; i < 100; ++i) {
      const auto j = std::uniform_int_distribution<std::uint64_t>{0, N - 1}(rand_eng);

      CHECK(*(first + j) == buff[j]);
      CHECK(*(last - j) == buff[N - j]);
    }
  }

  SECTION("subsequent calls to operator+=") {
    for (std::size_t i = 0; i < 10; ++i) {
      auto first = dset.begin<std::uint64_t>(32'000);
      auto last = dset.end<std::uint64_t>(32'000);
      std::size_t j = 0;

      while (first < last) {
        CHECK(*first == buff[j]);

        const auto step = std::uniform_int_distribution<std::size_t>{
            0, std::min(std::size_t(500), buff.size() - j)}(rand_eng);
        first += step;
        j += step;
      }
    }
  }

  SECTION("subsequent calls to operator-=") {
    for (std::size_t i = 0; i < 10; ++i) {
      auto first = dset.end<std::uint64_t>(32'000) - 1;
      auto last = dset.begin<std::uint64_t>(32'000);
      std::size_t j = buff.size() - 1;

      while (first > last) {
        CHECK(*first == buff[j]);

        const auto step =
            std::uniform_int_distribution<std::size_t>{0, std::min(std::size_t(500), j)}(rand_eng);

        first -= step;
        j -= step;
      }
    }
  }
}

}  // namespace hictk::cooler::test::dataset
