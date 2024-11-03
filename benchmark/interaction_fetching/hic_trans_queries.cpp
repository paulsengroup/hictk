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

static constexpr std::pair<std::string_view, std::string_view> range_small{
    "chr2L:15,000,000-15,100,000", "chrX:10,200,000-10,300,000"};

static constexpr std::pair<std::string_view, std::string_view> range_medium{
    "chr2L:5,000,000-6,000,000", "chrX:5,000,000-6,000,000"};

static constexpr std::pair<std::string_view, std::string_view> range_large{
    "chr2L:15,000,000-20,000,000", "chrX:15,000,000-20,000,000"};

template <typename N, bool sorted>
static void run_benchmark(const std::filesystem::path& path, std::uint32_t resolution,
                          std::string_view range1, std::string_view range2,
                          const balancing::Method& normalization) {
  BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}; {}bp; {}; {}"), range1, range2, resolution,
                                 sorted ? "sorted" : "unsorted",
                                 std::is_integral_v<N> ? "int" : "fp"))
  (Catch::Benchmark::Chronometer meter) {
    const hic::File hf(path.string(), resolution);
    meter.measure([&hf, &range1, &range2, &normalization]() {
      if constexpr (sorted) {
        return count_nnz<N>(hf, range1, range2, normalization);
      } else {
        return count_nnz_unsorted<N>(hf, range1, range2, normalization);
      }
    });
  };
}

TEST_CASE("hic::File::fetch (trans)") {
  const auto chroms = hic::File(test_file1.string(), resolutions.back()).chromosomes();

  for (const auto& path : {test_file1, test_file2}) {
    for (const auto& res : resolutions) {
      for (const auto& query : {range_small, range_medium, range_large}) {
        run_benchmark<std::uint32_t, true>(path, res, query.first, query.second,
                                           balancing::Method::NONE());
        run_benchmark<std::uint32_t, false>(path, res, query.first, query.second,
                                            balancing::Method::NONE());
        run_benchmark<double, true>(path, res, query.first, query.second, balancing::Method::KR());
        run_benchmark<double, false>(path, res, query.first, query.second, balancing::Method::KR());
      }
    }
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
