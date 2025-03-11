// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <string_view>

#include "./common.hpp"
#include "hictk/balancing/methods.hpp"
#include "hictk/benchmark/benchmark_installers.hpp"
#include "hictk/hic.hpp"

namespace hictk::benchmark {

// NOLINTBEGIN(*-avoid-magic-numbers, cert-err58-cpp)
static const TestCaseGenerator test_generator{
    "hic::File::fetch (gw)",
    std::array<std::string_view, 2>{"test/data/hic/4DNFIZ1ZVXC8.hic8",
                                    "test/data/hic/4DNFIZ1ZVXC8.hic9"},
    std::array<std::uint32_t, 4>{1000, 10000, 100000, 1000000},
    std::array<std::string_view, 1>{"GW"},
    std::array<std::string_view, 1>{"GW"},
    std::array<balancing::Method, 2>{balancing::Method::NONE(), balancing::Method::VC()}};
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp)

template <std::size_t I>
static void run_benchmark() {
  BENCHMARK_ADVANCED("benchmark")
  (Catch::Benchmark::Chronometer meter) {
    const auto& params = test_generator[I];
    const hic::File f(params.path.string(), params.resolution);
    if (params.normalization == balancing::Method::NONE()) {
      meter.measure([&]() {
        return count_nnz<std::uint32_t>(f, params.range1, params.range1, params.normalization);
      });
    } else {
      meter.measure([&]() {
        return count_nnz<double>(f, params.range1, params.range1, params.normalization);
      });
    }
  };
}

HICTK_REGISTER_BENCHMARKS(test_generator, run_benchmark)

void register_hic_gw_queries_benchmarks() { register_benchmarks(); }

}  // namespace hictk::benchmark
