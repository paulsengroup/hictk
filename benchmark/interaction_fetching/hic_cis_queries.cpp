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
#include "hictk/hic.hpp"

using namespace hictk;

// NOLINTBEGIN(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
static const std::filesystem::path test_file1{"test/data/hic/4DNFIZ1ZVXC8.hic8"};
static const std::filesystem::path test_file2{"test/data/hic/4DNFIZ1ZVXC8.hic9"};
static const std::vector<std::uint32_t> resolutions{1000,   5000,   10000,  25000,   50000,
                                                    100000, 250000, 500000, 1000000, 2500000};

static constexpr std::string_view range_small{"chr2L:5,000,000-5,100,000"};
static constexpr std::string_view range_medium{"chr2L:6,000,000-7,000,000"};
static constexpr std::string_view range_large{"chr2L:10,000,000-15,000,000"};

template <typename N, bool sorted>
static void run_benchmark(const std::filesystem::path& path, std::uint32_t resolution,
                          std::string_view range, const balancing::Method& normalization) {
  BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp; {}; {}"), range, resolution,
                                 sorted ? "sorted" : "unsorted",
                                 std::is_integral_v<N> ? "int" : "fp"))
  (Catch::Benchmark::Chronometer meter) {
    const hic::File hf(path.string(), resolution);
    meter.measure([&hf, &range, &normalization]() {
      if constexpr (sorted) {
        return count_nnz<N>(hf, range, range, normalization);
      } else {
        return count_nnz_unsorted<N>(hf, range, range, normalization);
      }
    });
  };
}

TEST_CASE("hic::File::fetch (cis)") {
  const auto chroms = hic::File(test_file1.string(), resolutions.back()).chromosomes();

  for (const auto& path : {test_file1, test_file2}) {
    for (const auto& res : resolutions) {
      for (const auto& range : {range_small, range_medium, range_large}) {
        run_benchmark<std::uint32_t, true>(path, res, range, balancing::Method::NONE());
        run_benchmark<std::uint32_t, false>(path, res, range, balancing::Method::NONE());
        run_benchmark<double, true>(path, res, range, balancing::Method::KR());
        run_benchmark<double, false>(path, res, range, balancing::Method::KR());
      }
    }
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
