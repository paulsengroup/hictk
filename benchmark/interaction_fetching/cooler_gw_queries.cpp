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

using namespace hictk;

// NOLINTBEGIN(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
static const std::filesystem::path test_file{"test/data/integration_tests/4DNFIZ1ZVXC8.mcool"};
static const std::vector<std::uint32_t> resolutions{1000,   5000,   10000,  25000,   50000,
                                                    100000, 250000, 500000, 1000000, 2500000};

TEST_CASE("cooler::File::fetch (gw; uint32)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file.string(), resolutions.back()))
                          .chromosomes();
  for (const auto& res : resolutions) {
    BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}bp"), res))
    (Catch::Benchmark::Chronometer meter) {
      const cooler::File clr(
          fmt::format(FMT_STRING("{}::/resolutions/{}"), test_file.string(), res));
      meter.measure([&clr]() {
        return count_nnz<std::uint32_t>(clr, 10'000'000, balancing::Method::NONE());
      });
    };
  }
}

TEST_CASE("cooler::File::fetch (gw; double)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file.string(), resolutions.back()))
                          .chromosomes();
  for (const auto& res : resolutions) {
    BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}bp"), res))
    (Catch::Benchmark::Chronometer meter) {
      const cooler::File clr(
          fmt::format(FMT_STRING("{}::/resolutions/{}"), test_file.string(), res));
      meter.measure(
          [&clr]() { return count_nnz<std::uint32_t>(clr, 10'000'000, balancing::Method::KR()); });
    };
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
