// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "./common.hpp"
#include "hictk/balancing/methods.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/file.hpp"

using namespace hictk;

// NOLINTBEGIN(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
static const std::filesystem::path test_file1{"test/data/integration_tests/4DNFIZ1ZVXC8.mcool"};
static const std::filesystem::path test_file2{"test/data/hic/4DNFIZ1ZVXC8.hic8"};
static const std::filesystem::path test_file3{"test/data/hic/4DNFIZ1ZVXC8.hic9"};
static const std::vector<std::uint32_t> resolutions{1000,   5000,   10000,  25000,   50000,
                                                    100000, 250000, 500000, 1000000, 2500000};

TEST_CASE("File::fetch (gw; uint32)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file1.string(), resolutions.back()))
                          .chromosomes();

  for (const auto& path : {test_file1, test_file2, test_file3}) {
    for (const auto& res : resolutions) {
      BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp"), path.extension(), res))
      (Catch::Benchmark::Chronometer meter) {
        const File f(path.string(), res);
        meter.measure(
            [&f]() { return count_nnz<std::uint32_t>(f, 10'000'000, balancing::Method::NONE()); });
      };
    }
  }
}

TEST_CASE("File::fetch (gw; double)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file1.string(), resolutions.back()))
                          .chromosomes();

  for (const auto& path : {test_file1, test_file2, test_file3}) {
    for (const auto& res : resolutions) {
      BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp"), path.extension(), res))
      (Catch::Benchmark::Chronometer meter) {
        const File f(path.string(), res);
        meter.measure(
            [&f]() { return count_nnz<std::uint32_t>(f, 10'000'000, balancing::Method::KR()); });
      };
    }
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
