

// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <random>
#include <string>
#include <string_view>
#include <utility>

#include "hictk/balancing/methods.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/reference.hpp"

struct Params {
  std::string_view label{};
  bool cis{};

  double avg_height{1.0e6};
  double avg_width{1.0e6};
  double height_std{250.0e3};
  double width_std{250.0e3};
  std::size_t num_queries{1};
  hictk::balancing::Method normalization{hictk::balancing::Method::NONE()};
  std::uint64_t seed{123456789};
};

[[nodiscard]] inline std::pair<std::string, std::string> generate_query(
    std::mt19937_64& rand_eng, const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2,
    double avg_height, double avg_width, double height_std, double width_std) {
  assert(chrom1 <= chrom2);

  const auto pos1 = std::uniform_int_distribution<std::uint32_t>{0U, chrom1.size() - 1}(rand_eng);
  const auto pos2 = std::uniform_int_distribution<std::uint32_t>{0U, chrom2.size() - 1}(rand_eng);

  const auto height = static_cast<std::uint32_t>(
      std::clamp(std::normal_distribution<double>{avg_height, height_std}(rand_eng), 1.0,
                 static_cast<double>(chrom1.size())));
  const auto width = static_cast<std::uint32_t>(
      std::clamp(std::normal_distribution<double>{avg_width, width_std}(rand_eng), 1.0,
                 static_cast<double>(chrom2.size())));

  auto start1 = height >= pos1 ? std::uint32_t{0} : pos1 - height;
  auto start2 = width >= pos2 ? std::uint32_t{0} : pos2 - width;

  if (chrom1 == chrom2 && start1 > start2) {
    std::swap(start1, start2);
  }

  auto end1 = std::min(start1 + height, chrom1.size());
  auto end2 = std::min(start2 + width, chrom2.size());

  return std::make_pair(fmt::format(FMT_STRING("{}:{}-{}"), chrom1.name(), start1, end1),
                        fmt::format(FMT_STRING("{}:{}-{}"), chrom2.name(), start2, end2));
}

[[nodiscard]] inline std::vector<std::pair<std::string, std::string>> generate_queries(
    const hictk::Chromosome& chrom1, const hictk::Chromosome& chrom2, std::size_t num_queries,
    double avg_height, double avg_width, double height_std, double width_std, std::uint64_t seed) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(seed == 0 ? rd() : seed);

  std::vector<std::pair<std::string, std::string>> queries(num_queries);

  std::generate(queries.begin(), queries.end(), [&]() {
    return generate_query(rand_eng, chrom1, chrom2, avg_height, avg_width, height_std, width_std);
  });

  return queries;
}

[[nodiscard]] inline std::discrete_distribution<std::uint32_t> init_chromosome_selector(
    const hictk::Reference& chroms) {
  std::vector<std::uint32_t> weights{};
  for (const auto& chrom : chroms) {
    if (chrom.is_all()) {
      weights.push_back(0);
      continue;
    }
    weights.push_back(chrom.size());
  }

  return {weights.begin(), weights.end()};
}

[[nodiscard]] inline std::vector<std::pair<std::string, std::string>> generate_queries_cis(
    const hictk::Reference& chroms, std::size_t num_queries, double avg_height, double avg_width,
    double height_std, double width_std, std::uint64_t seed) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(seed == 0 ? rd() : seed);

  auto chrom_selector = init_chromosome_selector(chroms);

  std::vector<std::pair<std::string, std::string>> queries(num_queries);

  std::generate(queries.begin(), queries.end(), [&]() {
    const auto& chrom = chroms.at(chrom_selector(rand_eng));
    return generate_query(rand_eng, chrom, chrom, avg_height, avg_width, height_std, width_std);
  });

  return queries;
}

[[nodiscard]] inline std::vector<std::pair<std::string, std::string>> generate_queries_trans(
    const hictk::Reference& chroms, std::size_t num_queries, double avg_height, double avg_width,
    double height_std, double width_std, std::uint64_t seed) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(seed == 0 ? rd() : seed);

  auto chrom_selector = init_chromosome_selector(chroms);

  std::vector<std::pair<std::string, std::string>> queries(num_queries);

  std::generate(queries.begin(), queries.end(), [&]() {
    const auto& chrom1 = chroms.at(chrom_selector(rand_eng));
    while (true) {
      const auto& chrom2 = chroms.at(chrom_selector(rand_eng));
      if (chrom1 == chrom2) {
        continue;
      }
      return generate_query(rand_eng, chrom1, chrom2, avg_height, avg_width, height_std, width_std);
    }
  });

  return queries;
}

template <typename N, typename File>
[[nodiscard]] inline std::ptrdiff_t count_nnz(const File& file, std::string_view range1,
                                              std::string_view range2,
                                              const hictk::balancing::Method& normalization) {
  const auto sel = file.fetch(range1, range2, normalization);

  return std::distance(sel.template begin<N>(), sel.template end<N>());
}
