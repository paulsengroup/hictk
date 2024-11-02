// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/cooler.hpp"

#include <fmt/format.h>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "./common.hpp"
#include "hictk/balancing/methods.hpp"

using namespace hictk;

static const std::filesystem::path test_file{"test/data/integration_tests/4DNFIZ1ZVXC8.mcool"};
static const std::vector<std::uint32_t> resolutions{1000,   5000,   10000,  25000,   50000,
                                                    100000, 250000, 500000, 1000000, 2500000};

static auto w = balancing::Method::KR();

static const std::vector<Params> params_uint{
    {"cis; small; normalization=NONE; symmetric", true, 100e3, 100e3, 25e3, 25e3},
    {"cis; medium; normalization=NONE; symmetric", true},
    {"cis; large; normalization=NONE; symmetric", true, 5e6, 5e6, 500e3, 500e3},
    {"trans; small; normalization=NONE; symmetric", false, 100e3, 100e3, 25e3, 25e3},
    {"trans; medium; normalization=NONE; symmetric", false},
    {"trans; large; normalization=NONE; symmetric", false, 5e6, 5e6, 500e3, 500e3},
};

static const std::vector<Params> params_fp{
    {"cis; small; normalization=weight; symmetric", true, 100e3, 100e3, 25e3, 25e3, 1, w},
    {"cis; medium; normalization=weight; symmetric", true, 1.0e6, 1.0e6, 250e3, 250e3, 1, w},
    {"cis; large; normalization=weight; symmetric", true, 5e6, 5e6, 500e3, 500e3, 1, w},
    {"trans; small; normalization=weight; symmetric", false, 100e3, 100e3, 25e3, 25e3, 1, w},
    {"trans; medium; normalization=weight; symmetric", false, 1.0e6, 1.0e6, 250e3, 250e3, 1, w},
    {"trans; large; normalization=weight; symmetric", false, 5e6, 5e6, 500e3, 500e3, 1, w},
};

TEST_CASE("cooler::File::fetch (uint32)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file.string(), resolutions.back()))
                          .chromosomes();
  for (const auto& res : resolutions) {
    for (const auto& [label, cis, avg_height, avg_width, height_std, width_std, num_queries,
                      normalization_, seed] : params_uint) {
      const auto normalization = normalization_;
      const auto& chrom1 = chroms.at(0);
      const auto& chrom2 = cis ? chrom1 : chroms.at(3);
      const auto queries = generate_queries(chrom1, chrom2, num_queries, avg_height, avg_width,
                                            height_std, width_std, seed);

      BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp"), label, res))
      (Catch::Benchmark::Chronometer meter) {
        const cooler::File clr(
            fmt::format(FMT_STRING("{}::/resolutions/{}"), test_file.string(), res));
        meter.measure([&clr, &queries, &normalization]() {
          std::int64_t nnz{};
          for (const auto& [range1, range2] : queries) {
            nnz += count_nnz<std::uint32_t>(clr, range1, range2, normalization);
          }
          return nnz;
        });
      };
    }
  }
}

TEST_CASE("cooler::File::fetch (double)") {
  const auto chroms = cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"),
                                               test_file.string(), resolutions.back()))
                          .chromosomes();
  for (const auto& res : resolutions) {
    for (const auto& [label, cis, avg_height, avg_width, height_std, width_std, num_queries,
                      normalization_, seed] : params_fp) {
      const auto normalization = normalization_;
      const auto& chrom1 = chroms.at(0);
      const auto& chrom2 = cis ? chrom1 : chroms.at(3);
      const auto queries = generate_queries(chrom1, chrom2, num_queries, avg_height, avg_width,
                                            height_std, width_std, seed);

      BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{}; {}bp"), label, res))
      (Catch::Benchmark::Chronometer meter) {
        const cooler::File clr(
            fmt::format(FMT_STRING("{}::/resolutions/{}"), test_file.string(), res));
        meter.measure([&clr, &queries, &normalization]() {
          std::int64_t nnz{};
          for (const auto& [range1, range2] : queries) {
            nnz += count_nnz<double>(clr, range1, range2, normalization);
          }
          return nnz;
        });
      };
    }
  }
}
