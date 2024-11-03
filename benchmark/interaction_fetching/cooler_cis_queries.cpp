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

template <typename N>
static void run_benchmark(const std::filesystem::path& path, std::uint32_t resolution,
                          std::string_view range, const balancing::Method& normalization) {
  BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp; {}"), range, resolution,
                                 std::is_integral_v<N> ? "int" : "fp"))
  (Catch::Benchmark::Chronometer meter) {
    const cooler::File clr(
        fmt::format(FMT_STRING("{}::/resolutions/{}"), path.string(), resolution));
    meter.measure([&clr, &range, &normalization]() {
      return count_nnz<N>(clr, range, range, normalization);
    });
  };
}

TEST_CASE("cooler::File::fetch (cis)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file.string(), resolutions.back()))
                          .chromosomes();
  for (const auto& res : resolutions) {
    for (const auto& range : {range_small, range_medium, range_large}) {
      run_benchmark<std::uint32_t>(test_file, res, range, balancing::Method::NONE());
      run_benchmark<double>(test_file, res, range, balancing::Method::KR());
    }
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
