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
#include "hictk/hic.hpp"

using namespace hictk;

// NOLINTBEGIN(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
static const std::filesystem::path test_file1{"test/data/hic/4DNFIZ1ZVXC8.hic8"};
static const std::filesystem::path test_file2{"test/data/hic/4DNFIZ1ZVXC8.hic9"};
static const std::vector<std::uint32_t> resolutions{1000,   5000,   10000,  25000,   50000,
                                                    100000, 250000, 500000, 1000000, 2500000};

static const auto w = balancing::Method::KR();

static const std::vector<Params> params_uint{
    {"trans; small; normalization=NONE; symmetric", false, 100e3, 100e3, 25e3, 25e3},
    {"trans; medium; normalization=NONE; symmetric", false},
    {"trans; large; normalization=NONE; symmetric", false, 5e6, 5e6, 500e3, 500e3},
};

static const std::vector<Params> params_fp{
    {"trans; small; normalization=weight; symmetric", false, 100e3, 100e3, 25e3, 25e3, 1, w},
    {"trans; medium; normalization=weight; symmetric", false, 1.0e6, 1.0e6, 250e3, 250e3, 1, w},
    {"trans; large; normalization=weight; symmetric", false, 5e6, 5e6, 500e3, 500e3, 1, w},
};

TEST_CASE("hic::File::fetch (trans; uint32)") {
  const auto chroms = hic::File(test_file1.string(), resolutions.back()).chromosomes();

  for (const auto& path : {test_file1, test_file2}) {
    for (const auto& res : resolutions) {
      for (const auto& [label, cis, avg_height, avg_width, height_std, width_std, num_queries,
                        normalization_, seed] : params_uint) {
        const auto normalization = normalization_;
        const auto& chrom1 = chroms.at(1);
        const auto& chrom2 = cis ? chrom1 : chroms.at(4);
        const auto queries = generate_queries(chrom1, chrom2, num_queries, avg_height, avg_width,
                                              height_std, width_std, seed);

        BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp"), label, res))
        (Catch::Benchmark::Chronometer meter) {
          const hic::File hf(path.string(), res);
          meter.measure([&hf, &queries, &normalization]() {
            std::int64_t nnz{};
            for (const auto& [range1, range2] : queries) {
              nnz += count_nnz<std::uint32_t>(hf, range1, range2, normalization);
            }
            return nnz;
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
      for (const auto& [label, cis, avg_height, avg_width, height_std, width_std, num_queries,
                        normalization_, seed] : params_fp) {
        const auto normalization = normalization_;
        const auto& chrom1 = chroms.at(1);
        const auto& chrom2 = cis ? chrom1 : chroms.at(4);
        const auto queries = generate_queries(chrom1, chrom2, num_queries, avg_height, avg_width,
                                              height_std, width_std, seed);

        BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp"), label, res))
        (Catch::Benchmark::Chronometer meter) {
          const hic::File hf(path.string(), res);
          meter.measure([&hf, &queries, &normalization]() {
            std::int64_t nnz{};
            for (const auto& [range1, range2] : queries) {
              nnz += count_nnz<double>(hf, range1, range2, normalization);
            }
            return nnz;
          });
        };
      }
    }
  }
}
// NOLINTEND(*-avoid-magic-numbers, cert-err58-cpp, readability-function-cognitive-complexity)
