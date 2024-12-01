// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/benchmark/catch_constructor.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <iterator>
#include <random>
#include <vector>

#include "hictk/benchmark/hg38.hpp"
#include "hictk/bin.hpp"
#include "hictk/bin_table_fixed.hpp"
#include "hictk/chromosome.hpp"

namespace hictk::benchmark {

// NOLINTBEGIN(*-avoid-magic-numbers)
[[nodiscard]] static std::vector<Bin> generate_bins(std::size_t size, bool erase_ids) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  const BinTableFixed bin_table({hg38.begin(), hg38.end()}, 1'000);
  std::vector<Bin> bins;
  bins.reserve(size);

  while (bins.size() != size) {
    std::sample(bin_table.begin(), bin_table.end(), std::back_inserter(bins), size - bins.size(),
                rand_eng);
  }

  if (erase_ids) {
    std::transform(bins.begin(), bins.end(), bins.begin(),
                   [](const auto& bin) { return Bin{bin.chrom(), bin.start(), bin.end()}; });
  }

  return bins;
}

TEST_CASE("Bin") {
  BENCHMARK_ADVANCED("Construction")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = static_cast<std::size_t>(meter.runs());
    std::vector<Catch::Benchmark::storage_for<Bin>> storage(num_runs);

    const Chromosome chrom{0, "chr1", 123'456'789};

    meter.measure(
        [&storage, &chrom](std::size_t i) { storage[i].construct(chrom, 10'000'000, 11'000'000); });
  };

  BENCHMARK_ADVANCED("Destruction")
  (Catch::Benchmark::Chronometer meter) {
    const auto num_runs = static_cast<std::size_t>(meter.runs());
    std::vector<Catch::Benchmark::destructable_object<Bin>> storage(num_runs);

    const Chromosome chrom{0, "chr1", 123'456'789};

    for (auto& bin : storage) {
      bin.construct(chrom, 10'000'000, 11'000'000);
    }

    meter.measure([&storage](std::size_t i) { storage[i].destruct(); });
  };

  BENCHMARK_ADVANCED("sorting w/ id")
  (Catch::Benchmark::Chronometer meter) {
    const auto bins = generate_bins(1'000'000, false);
    std::vector<std::vector<Bin>> data(static_cast<std::size_t>(meter.runs()), bins);
    meter.measure([&data](std::size_t i) {
      std::size_t num_ops{};
      std::sort(data[i].begin(), data[i].end(), [&](const auto& bin1, const auto& bin2) {
        ++num_ops;
        return bin1 < bin2;
      });
      return num_ops;
    });
  };

  BENCHMARK_ADVANCED("sorting wo/ id")
  (Catch::Benchmark::Chronometer meter) {
    const auto bins = generate_bins(1'000'000, true);
    std::vector<std::vector<Bin>> data(static_cast<std::size_t>(meter.runs()), bins);
    meter.measure([&data](std::size_t i) {
      std::size_t num_ops{};
      std::sort(data[i].begin(), data[i].end(), [&](const auto& bin1, const auto& bin2) {
        ++num_ops;
        return bin1 < bin2;
      });
      return num_ops;
    });
  };
}
// NOLINTEND(*-avoid-magic-numbers)

}  // namespace hictk::benchmark
