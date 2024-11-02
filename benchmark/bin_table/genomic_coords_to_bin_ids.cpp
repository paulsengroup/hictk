// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <hictk/bin_table.hpp>
#include <hictk/bin_table_fixed.hpp>
#include <vector>

#include "./common.hpp"

using namespace hictk;

TEST_CASE("BinTable::at(chrom, pos)") {
  const std::vector<std::uint32_t> resolutions{10, 100, 1'000, 10'000, 100'000, 1'000'000};

  for (const auto &res : resolutions) {
    BENCHMARK_ADVANCED(fmt::format(FMT_STRING("hg38 ({}bp)"), res))
    (Catch::Benchmark::Chronometer meter) {
      const BinTable bin_table{hg38.begin(), hg38.end(), res};
      const auto coords =
          generate_genomic_coords(bin_table, static_cast<std::size_t>(meter.runs()));

      meter.measure([&bin_table, &coords](std::size_t i) {
        return bin_table.at(coords[i].first, coords[i].second);
      });
    };
  }
}

TEST_CASE("BinTableFixed::at(chrom, pos)") {
  const std::vector<std::uint32_t> resolutions{10, 100, 1'000, 10'000, 100'000, 1'000'000};

  for (const auto &res : resolutions) {
    BENCHMARK_ADVANCED(fmt::format(FMT_STRING("hg38 ({}bp)"), res))
    (Catch::Benchmark::Chronometer meter) {
      const BinTableFixed bin_table{hg38.begin(), hg38.end(), res};
      const auto coords =
          generate_genomic_coords(bin_table, static_cast<std::size_t>(meter.runs()));

      meter.measure([&bin_table, &coords](std::size_t i) {
        return bin_table.at(coords[i].first, coords[i].second);
      });
    };
  }
}

TEST_CASE("BinTableVariable::at(chrom, pos)") {
  const std::vector<std::uint32_t> resolutions{5'000, 10'000, 100'000, 1'000'000};

  for (const auto &res : resolutions) {
    BENCHMARK_ADVANCED(fmt::format(FMT_STRING("hg38 ({}bp)"), res))
    (Catch::Benchmark::Chronometer meter) {
      const auto bin_table = generate_variable_bin_table(res);
      const auto coords =
          generate_genomic_coords(bin_table, static_cast<std::size_t>(meter.runs()));

      meter.measure([&bin_table, &coords](std::size_t i) {
        return bin_table.at(coords[i].first, coords[i].second);
      });
    };
  }
}
