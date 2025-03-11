// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <string_view>
#include <vector>

#include "hictk/benchmark/hg38.hpp"
#include "hictk/bin_table_variable.hpp"
#include "hictk/chromosome.hpp"

namespace hictk::benchmark {

template <typename BinTable>
[[nodiscard]] inline std::vector<std::uint64_t> generate_bin_ids(const BinTable &bins,
                                                                 std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  std::vector<std::uint64_t> buff(size);
  std::generate(buff.begin(), buff.end(), [&]() {
    return std::uniform_int_distribution<std::uint64_t>{0, bins.size() - 1}(rand_eng);
  });

  return buff;
}

template <typename BinTable>
[[nodiscard]] inline auto generate_genomic_coords(const BinTable &bins, std::size_t size) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  using Coord = std::pair<std::uint32_t, std::uint32_t>;
  std::vector<Coord> buff(size);
  std::generate(buff.begin(), buff.end(), [&]() {
    const auto bin_id = std::uniform_int_distribution<std::uint64_t>{0, bins.size() - 1}(rand_eng);

    const auto chrom = bins.at(bin_id).chrom();
    const auto pos = std::uniform_int_distribution<std::uint32_t>{0, chrom.size() - 1}(rand_eng);

    return std::make_pair(chrom.id(), pos);
  });

  return buff;
}

[[nodiscard]] inline hictk::BinTableVariable<std::uint32_t> generate_variable_bin_table(
    std::uint32_t target_resolution) {
  std::random_device rd{};
  std::mt19937_64 rand_eng(rd());

  const auto resolution_avg = static_cast<double>(target_resolution);
  const auto resolution_std = std::max(10.0, resolution_avg / 10);

  auto generate_bin_size = [&](const hictk::Chromosome &chrom, std::uint32_t pos) {
    const auto bin_size =
        std::normal_distribution<double>{resolution_avg, resolution_std}(rand_eng);
    return static_cast<std::uint32_t>(
        std::clamp(bin_size, 1.0, static_cast<double>(chrom.size() - pos)));
  };

  std::vector<std::uint32_t> start_pos{};
  std::vector<std::uint32_t> end_pos{};

  for (const auto &chrom : hg38) {
    start_pos.push_back(0);
    end_pos.push_back(start_pos.back() + generate_bin_size(chrom, start_pos.back()));
    while (end_pos.back() < chrom.size()) {
      start_pos.push_back(end_pos.back());
      end_pos.push_back(start_pos.back() + generate_bin_size(chrom, start_pos.back()));
    }
  }

  return {hictk::Reference{hg38.begin(), hg38.end()}, start_pos, end_pos};
}

}  // namespace hictk::benchmark
