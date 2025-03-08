// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./init_bin_table.hpp"

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/reference.hpp"

namespace hictk::tools {

BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes, std::uint32_t bin_size) {
  auto chroms = Reference::from_chrom_sizes(path_to_chrom_sizes);
  return {chroms, bin_size};
}

struct ChromName {
  std::uint32_t id{};
  std::string name{};

  ChromName() = default;
  ChromName(std::uint32_t id_, std::string name_) : id(id_), name(std::move(name_)) {}

  [[nodiscard]] bool operator<(const ChromName& other) const { return id < other.id; }
};

struct ChromNameCmp {
  using is_transparent = int;

  [[nodiscard]] bool operator()(const std::shared_ptr<const ChromName>& chrom1,
                                const std::string& chrom2_name) const noexcept {
    return chrom1->name < chrom2_name;
  }
  [[nodiscard]] bool operator()(const std::string& chrom1_name,
                                const std::shared_ptr<const ChromName>& chrom2) const noexcept {
    return chrom1_name < chrom2->name;
  }
  [[nodiscard]] bool operator()(const std::string& chrom1_name,
                                const std::string& chrom2_name) const noexcept {
    return chrom1_name < chrom2_name;
  }
  [[nodiscard]] bool operator()(const std::shared_ptr<const ChromName>& chrom1,
                                const std::shared_ptr<const ChromName>& chrom2) const noexcept {
    return chrom1->name < chrom2->name;
  }
};

BinTable init_bin_table(const std::filesystem::path& path_to_bin_table) {
  std::ifstream ifs{};
  ifs.exceptions(std::ios::badbit);
  ifs.open(path_to_bin_table);

  phmap::btree_map<std::shared_ptr<const ChromName>, std::uint32_t, ChromNameCmp> chrom_sizes;
  std::vector<std::shared_ptr<const ChromName>> chrom_names{};
  std::vector<std::uint32_t> start_pos{};
  std::vector<std::uint32_t> end_pos{};

  std::string line{};
  std::uint32_t bin_size = 0;
  while (std::getline(ifs, line)) {
    auto [chrom, start, end] = GenomicInterval::parse_bed(line);
    auto match = chrom_sizes.find(chrom);
    if (match == chrom_sizes.end()) {
      const auto chrom_id = static_cast<std::uint32_t>(chrom_sizes.size());
      auto [it, _] =
          chrom_sizes.emplace(std::make_shared<const ChromName>(chrom_id, std::move(chrom)), 0);
      match = it;
    }

    chrom_names.emplace_back(match->first);
    start_pos.push_back(start);
    end_pos.push_back(end);

    bin_size = std::max(bin_size, end - start);

    match->second = std::max(match->second, end);
  }

  if (chrom_sizes.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to import bins from \"{}\": file appears to be empty"),
                    path_to_bin_table));
  }

  assert(bin_size != 0);

  bool fixed_bin_size = true;
  for (std::size_t i = 0; i < chrom_names.size(); ++i) {
    const auto record_span = end_pos[i] - start_pos[i];
    if (record_span != bin_size) {
      const auto chrom_size = chrom_sizes.find(chrom_names[i])->second;
      if (record_span != chrom_size) {
        fixed_bin_size = false;
        break;
      }
    }
  }

  std::vector<Chromosome> chroms(chrom_sizes.size());
  std::transform(chrom_sizes.begin(), chrom_sizes.end(), chroms.begin(),
                 [](const auto& kv) -> Chromosome {
                   const auto& [chrom_id, chrom_name] = *kv.first;
                   const auto chrom_size = kv.second;
                   return {chrom_id, chrom_name, chrom_size};
                 });
  std::sort(chroms.begin(), chroms.end());

  if (fixed_bin_size) {
    SPDLOG_INFO("detected bin table with uniform bin size.");
    return {chroms.begin(), chroms.end(), bin_size};
  }

  SPDLOG_INFO("detected bin table with variable bin size.");
  return {Reference{chroms.begin(), chroms.end()}, start_pos, end_pos};
}

BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes,
                        const std::filesystem::path& path_to_bin_table, std::uint32_t bin_size) {
  if (!path_to_bin_table.empty()) {
    return init_bin_table(path_to_bin_table);
  }
  return init_bin_table(path_to_chrom_sizes, bin_size);
}

}  // namespace hictk::tools
