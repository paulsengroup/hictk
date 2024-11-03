// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <string_view>
#include <utility>
#include <vector>

#include "./common.hpp"
#include "hictk/balancing/methods.hpp"
#include "hictk/cooler/cooler.hpp"

using namespace hictk;

// NOLINTBEGIN(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
static const std::filesystem::path test_file{"test/data/integration_tests/4DNFIZ1ZVXC8.mcool"};
static const std::vector<std::uint32_t> resolutions{1000,   5000,   10000,  25000,   50000,
                                                    100000, 250000, 500000, 1000000, 2500000};

static constexpr std::string_view range_small{"chr2L:5,000,000-5,100,000"};
static constexpr std::string_view range_medium{"chr2L:6,000,000-7,000,000"};
static constexpr std::string_view range_large{"chr2L:10,000,000-15,000,000"};

TEST_CASE("cooler::File::fetch (cis; uint32)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file.string(), resolutions.back()))
                          .chromosomes();
  for (const auto& res : resolutions) {
    for (const auto& range : {range_small, range_medium, range_large}) {
      BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp"), range, res))
      (Catch::Benchmark::Chronometer meter) {
        const cooler::File clr(
            fmt::format(FMT_STRING("{}::/resolutions/{}"), test_file.string(), res));
        meter.measure([&clr, &range]() {
          return count_nnz<std::uint32_t>(clr, range, range, balancing::Method::NONE());
        });
      };
    }
  }
}

TEST_CASE("cooler::File::fetch (cis; double)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file.string(), resolutions.back()))
                          .chromosomes();
  for (const auto& res : resolutions) {
    for (const auto& range : {range_small, range_medium, range_large}) {
      BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp"), range, res))
      (Catch::Benchmark::Chronometer meter) {
        const cooler::File clr(
            fmt::format(FMT_STRING("{}::/resolutions/{}"), test_file.string(), res));
        meter.measure([&clr, &range]() {
          return count_nnz<double>(clr, range, range, balancing::Method::KR());
        });
      };
    }
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
