// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>
#include <random>
#include <set>

#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::dataset {
TEST_CASE("Cooler: dataset large read/write", "[dataset][long]") {
  const auto path = testdir() / "test_dataset_large_rw.h5";

  constexpr std::uint64_t seed{4195331987557451569};
  constexpr std::size_t N = 5'000'000;
  {
    const RootGroup grp{HighFive::File(path.string(), HighFive::File::Truncate).getGroup("/")};
    Dataset dset{grp, "int", std::uint8_t{}};

    std::vector<std::uint8_t> buff(1'000'000);
    std::mt19937_64 rand_eng{seed};  // NOLINT(cert-msc32-c,cert-msc51-cpp)
    while (dset.size() != N) {
      std::generate(buff.begin(), buff.end(),
                    [&]() { return static_cast<std::uint8_t>(rand_eng()); });
      dset.append(buff);
    }
    CHECK(dset.size() == N);
  }

  const RootGroup grp{HighFive::File(path.string()).getGroup("/")};
  const Dataset dset{grp, "int"};
  REQUIRE(dset.size() == N);

  std::mt19937_64 rand_eng{seed};  // NOLINT(cert-msc32-c,cert-msc51-cpp)
  std::for_each(dset.begin<std::uint8_t>(256'000), dset.end<std::uint8_t>(256'000),
                [&](const auto& n) { CHECK(n == static_cast<std::uint8_t>(rand_eng())); });
}

}  // namespace hictk::cooler::test::dataset
