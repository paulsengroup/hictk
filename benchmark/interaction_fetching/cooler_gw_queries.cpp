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

template <typename N>
static void run_benchmark(const std::filesystem::path& path, std::uint32_t resolution,
                          const balancing::Method& normalization) {
  BENCHMARK_ADVANCED(
      fmt::format(FMT_STRING("{}bp; {}"), resolution, std::is_integral_v<N> ? "int" : "fp"))
  (Catch::Benchmark::Chronometer meter) {
    const cooler::File clr(
        fmt::format(FMT_STRING("{}::/resolutions/{}"), path.string(), resolution));
    meter.measure(
        [&clr, &normalization]() { return count_nnz<N>(clr, 10'000'000, normalization); });
  };
}

TEST_CASE("cooler::File::fetch (gw)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file.string(), resolutions.back()))
                          .chromosomes();
  for (const auto& res : resolutions) {
    run_benchmark<std::uint32_t>(test_file, res, balancing::Method::NONE());
    run_benchmark<double>(test_file, res, balancing::Method::KR());
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
