// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <random>
#include <string_view>
#include <vector>

#include "hictk/bin_table_variable.hpp"
#include "hictk/chromosome.hpp"

// clang-format off
// NOLINTNEXTLINE(cert-err58-cpp)
inline const std::vector hg38{
        hictk::Chromosome{0,  "chr1",  248956422},
        hictk::Chromosome{1,  "chr2",  242193529},
        hictk::Chromosome{2,  "chr3",  198295559},
        hictk::Chromosome{3,  "chr4",  190214555},
        hictk::Chromosome{4,  "chr5",  181538259},
        hictk::Chromosome{5,  "chr6",  170805979},
        hictk::Chromosome{6,  "chr7",  159345973},
        hictk::Chromosome{7,  "chr8",  145138636},
        hictk::Chromosome{8,  "chr9",  138394717},
        hictk::Chromosome{9,  "chr10", 133797422},
        hictk::Chromosome{10, "chr11", 135086622},
        hictk::Chromosome{11, "chr12", 133275309},
        hictk::Chromosome{12, "chr13", 114364328},
        hictk::Chromosome{13, "chr14", 107043718},
        hictk::Chromosome{14, "chr15", 101991189},
        hictk::Chromosome{15, "chr16", 90338345},
        hictk::Chromosome{16, "chr17", 83257441},
        hictk::Chromosome{17, "chr18", 80373285},
        hictk::Chromosome{18, "chr19", 58617616},
        hictk::Chromosome{19, "chr20", 64444167},
        hictk::Chromosome{20, "chr21", 46709983},
        hictk::Chromosome{21, "chr22", 50818468},
        hictk::Chromosome{22, "chrX",  156040895},
        hictk::Chromosome{23, "chrY",  57227415}
};
// clang-format on

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

  using Coord = std::pair<std::string_view, std::uint32_t>;
  std::vector<Coord> buff(size);
  std::generate(buff.begin(), buff.end(), [&]() {
    const auto bin_id = std::uniform_int_distribution<std::uint64_t>{0, bins.size() - 1}(rand_eng);

    const auto chrom = bins.at(bin_id).chrom();
    const auto pos = std::uniform_int_distribution<std::uint32_t>{0, chrom.size() - 1}(rand_eng);

    return std::make_pair(chrom.name(), pos);
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
