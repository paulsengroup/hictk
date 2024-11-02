// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cstddef>
#include <cstdint>
#include <random>
#include <vector>

#include "hictk/reference.hpp"

using namespace hictk;

[[nodiscard]] static Reference generate_reference(std::size_t num_chroms) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<std::string> names(num_chroms);
  std::vector<std::uint32_t> sizes(num_chroms);

  for (std::size_t i = 0; i < num_chroms; ++i) {
    names[i] = fmt::format(FMT_STRING("chr{}"), i + 1);
    sizes[i] = std::uniform_int_distribution<std::uint32_t>{1'000'000, 500'000'000}(rand_eng);
  }

  return {names.begin(), names.end(), sizes.begin()};
}

[[nodiscard]] static std::vector<std::string_view> generate_chrom_names(const Reference& chroms,
                                                                        std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<std::string_view> names(size);
  std::generate(names.begin(), names.end(), [&]() {
    const auto chrom_id = std::uniform_int_distribution<std::uint32_t>{
        0, static_cast<std::uint32_t>(chroms.size() - 1)}(rand_eng);
    return chroms.at(chrom_id).name();
  });

  return names;
}

[[nodiscard]] static std::vector<std::uint32_t> generate_chrom_ids(const Reference& chroms,
                                                                   std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<std::uint32_t> ids(size);
  std::generate(ids.begin(), ids.end(), [&]() {
    return std::uniform_int_distribution<std::uint32_t>{
        0, static_cast<std::uint32_t>(chroms.size() - 1)}(rand_eng);
  });

  return ids;
}

TEST_CASE("Reference::at(name)") {
  const std::vector<std::size_t> num_chroms{5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000};

  for (const auto& size : num_chroms) {
    BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{} chromosomes)"), size))
    (Catch::Benchmark::Chronometer meter) {
      const auto chroms = generate_reference(size);
      const auto names = generate_chrom_names(chroms, static_cast<std::size_t>(meter.runs()));

      meter.measure([&chroms, &names](std::size_t i) { return chroms.at(names[i]); });
    };
  }
}

TEST_CASE("Reference::at(id)") {
  const std::vector<std::size_t> num_chroms{5, 10, 20, 30, 40, 50, 100, 200, 300, 400, 500, 1000};

  for (const auto& size : num_chroms) {
    BENCHMARK_ADVANCED(fmt::format(FMT_STRING("{} chromosomes)"), size))
    (Catch::Benchmark::Chronometer meter) {
      const auto chroms = generate_reference(size);
      const auto ids = generate_chrom_ids(chroms, static_cast<std::size_t>(meter.runs()));

      meter.measure([&chroms, &ids](std::size_t i) { return chroms.at(ids[i]); });
    };
  }
}
