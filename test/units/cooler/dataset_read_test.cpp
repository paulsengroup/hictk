// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <numeric>
#include <random>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/test/testdir.hpp"
#include "hictk/variant_buff.hpp"

namespace hictk::cooler::test::dataset {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Cooler: dataset read", "[dataset][short]") {
  const auto path = datadir / "cooler_test_file.cool";
  const RootGroup grp{HighFive::File(path.string()).getGroup("/")};

  SECTION("fixed str") {
    SECTION("vector") {
      constexpr std::array<std::string_view, 3> expected{"1", "2", "3"};
      std::vector<std::string> buff{};

      Dataset{grp, "chroms/name"}.read(buff, expected.size());
      REQUIRE(buff.size() == 3);

      for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }
    }

    SECTION("atomic") {
      const Dataset dset{grp, "chroms/name"};
      std::string buff;
      dset.read(buff, 9);
      CHECK(buff == "10");
      CHECK(dset.read_last<std::string>() == "X");
      CHECK(std::get<std::string>(dset.read_last()) == "X");
    }
  }

  SECTION("numeric") {
    using T = std::int32_t;
    constexpr std::array<T, 10> expected{0,       100'000, 200'000, 300'000, 400'000,
                                         500'000, 600'000, 700'000, 800'000, 900'000};

    constexpr std::size_t nnz_expected = 107'041;
    constexpr std::int32_t sum_expected = 395'465;

    SECTION("vector<T>") {
      std::vector<T> buff{};
      std::ignore = Dataset{grp, "bins/start"}.read(buff, expected.size());

      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }

      Dataset{grp, "pixels/count"}.read_all(buff);
      CHECK(buff.size() == nnz_expected);
      CHECK(std::accumulate(buff.begin(), buff.end(), 0) == sum_expected);
    }

    SECTION("variant buff") {
      hictk::internal::VariantBuffer vbuff{std::size_t{0}, 0.0};
      std::ignore = Dataset{grp, "bins/start"}.read(vbuff, expected.size());
      const auto& buff = vbuff.get<T>();
      REQUIRE(buff.size() == expected.size());
      for (std::size_t i = 0; i < buff.size(); ++i) {
        CHECK(buff[i] == expected[i]);
      }

      Dataset{grp, "pixels/count"}.read_all(vbuff);
      CHECK(vbuff.size<T>() == nnz_expected);
      CHECK(std::accumulate(vbuff.begin<T>(), vbuff.end<T>(), 0) == sum_expected);
    }

    SECTION("atomic") {
      const Dataset dset{grp, "chroms/length"};
      std::uint64_t buff{};
      dset.read(buff, 2);
      CHECK(buff == 159'599'783);

      CHECK(dset.read_last<std::int32_t>() == 166'650'296);
      CHECK(std::get<std::int32_t>(dset.read_last()) == 166'650'296);
    }

    SECTION("enum") { CHECK(Dataset{grp, "bins/chrom"}.read<std::uint32_t>(0) == 0); }
  }
}

[[nodiscard]] static std::vector<std::int16_t> generate_sorted_vec_of_unique_numbers(
    std::mt19937_64& rand_eng) {
  // Generate a random vector of unique numbers sorted in ascending order
  std::vector<std::int16_t> buff(10'000);
  std::generate(buff.begin(), buff.end(), [&]() {
    return std::uniform_int_distribution<std::int16_t>{-10'000, 10'000}(rand_eng);
  });
  std::sort(buff.begin(), buff.end());
  buff.erase(std::unique(buff.begin(), buff.end()), buff.end());

  return buff;
}

TEST_CASE("Cooler: dataset lower_bound", "[dataset][medium]") {
  const auto path = testdir() / "test_dataset_lower_bound.h5";

  constexpr std::uint64_t seed{18125230607725213391ULL};
  std::mt19937_64 rand_eng{seed};  // NOLINT(cert-msc32-c,cert-msc51-cpp)

  const auto buff = generate_sorted_vec_of_unique_numbers(rand_eng);
  {
    // write vector to .hdf5 file
    const RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
    Dataset{grp, "data", std::int16_t{}}.append(buff);
  }

  const RootGroup grp{HighFive::File(path.string(), HighFive::File::ReadOnly).getGroup("/")};
  const Dataset dset{grp, "data"};

  constexpr std::ptrdiff_t chunk_size = 123;
  REQUIRE(buff.size() / chunk_size > 50);  // make sure we have a good number of chunks to play with

  SECTION("within (value present)") {
    const std::int16_t value = buff[buff.size() / 2];
    auto first = dset.begin<std::int16_t>(chunk_size);
    auto last = dset.end<std::int16_t>(0);

    auto found = Dataset::lower_bound(first, last, value);
    REQUIRE(found != last);
    CHECK(*found == value);
  }

  SECTION("within (value missing)") {
    std::size_t i = buff.size() / 2;
    REQUIRE(i > 0);
    for (; i < buff.size(); ++i) {
      const auto v0 = buff[i - 1];
      const auto v1 = buff[i];
      if (v1 - v0 > 1) {
        // found a gap in the sequence of numbers
        break;
      }
    }

    REQUIRE(i != buff.size());
    const std::int16_t value = buff[i - 1] + 1;  // NOLINT(*-narrowing-conversions)
    const auto next_value = buff[i];

    auto first = dset.begin<std::int16_t>(chunk_size);
    auto last = dset.end<std::int16_t>(0);

    auto found = Dataset::lower_bound(first, last, value);
    REQUIRE(found != last);
    CHECK(*found == next_value);
  }

  SECTION("upstream") {
    const std::int16_t value = buff.front() - 1;  // NOLINT(*-narrowing-conversions)
    auto first = dset.begin<std::int16_t>(chunk_size);
    auto last = dset.end<std::int16_t>(0);

    CHECK(Dataset::lower_bound(first, last, value) == first);
  }

  SECTION("downstream") {
    const std::int16_t value = buff.back() + 1;  // NOLINT(*-narrowing-conversions)
    auto first = dset.begin<std::int16_t>(chunk_size);
    auto last = dset.end<std::int16_t>(0);

    CHECK(Dataset::lower_bound(first, last, value) == last);
  }

  SECTION("randomized") {
    const auto max_offset1 = static_cast<std::ptrdiff_t>(buff.size() - (buff.size() / 2));
    const auto max_offset2 = static_cast<std::ptrdiff_t>(buff.size());

    for (std::size_t i = 0; i < 25'000; ++i) {
      const auto offset1 = std::uniform_int_distribution<std::ptrdiff_t>{0, max_offset1}(rand_eng);
      const auto offset2 = std::uniform_int_distribution{offset1, max_offset2}(rand_eng);

      const auto guess_inside = std::bernoulli_distribution{0.9}(rand_eng);
      const auto lb = guess_inside ? buff[static_cast<std::size_t>(offset1)]
                                   : std::numeric_limits<std::int16_t>::min();
      const auto ub = guess_inside ? buff[static_cast<std::size_t>(offset2)]
                                   : std::numeric_limits<std::int16_t>::max();
      const auto value = std::uniform_int_distribution{lb, std::max(lb, ub)}(rand_eng);

      const auto expected = std::lower_bound(buff.begin() + offset1, buff.begin() + offset2, value);
      const auto found = Dataset::lower_bound(dset.begin<std::int16_t>(chunk_size) + offset1,
                                              dset.begin<std::int16_t>(0) + offset2, value);

      if (expected == buff.end()) {
        CHECK(found == dset.end<std::int16_t>());
      } else if (found == dset.end<std::int16_t>()) {
        CHECK(expected == buff.end());
      } else {
        CHECK(*found == *expected);
      }
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::dataset
