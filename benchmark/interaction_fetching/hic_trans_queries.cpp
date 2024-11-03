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

TEST_CASE("hic::File::fetch (trans; uint32)") {
  const auto chroms = hic::File(test_file1.string(), resolutions.back()).chromosomes();

  for (const auto& path : {test_file1, test_file2}) {
    for (const auto& res : resolutions) {
      for (const auto& query : {range_small, range_medium, range_large}) {
        BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}; {}bp"), query.first, query.second, res))
        (Catch::Benchmark::Chronometer meter) {
          const hic::File hf(path.string(), res);
          meter.measure([&hf, range1 = query.first, range2 = query.second]() {
            return count_nnz<std::uint32_t>(hf, range1, range2, balancing::Method::NONE());
          });
        };
      }
    }
  }
}

TEST_CASE("hic::File::fetch (trans; double)") {
  const auto chroms = hic::File(test_file1.string(), resolutions.back()).chromosomes();

  for (const auto& path : {test_file1, test_file2}) {
    for (const auto& res : resolutions) {
      for (const auto& query : {range_small, range_medium, range_large}) {
        BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}; {}bp"), query.first, query.second, res))
        (Catch::Benchmark::Chronometer meter) {
          const hic::File hf(path.string(), res);
          meter.measure([&hf, range1 = query.first, range2 = query.second]() {
            return count_nnz<double>(hf, range1, range2, balancing::Method::KR());
          });
        };
      }
    }
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
