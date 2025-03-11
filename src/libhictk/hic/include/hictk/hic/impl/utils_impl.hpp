// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/file_reader.hpp"

namespace hictk::hic::utils {
inline std::vector<std::uint32_t> list_resolutions(const std::filesystem::path& path, bool sorted) {
  auto resolutions = hic::internal::HiCFileReader(path.string()).header().resolutions;
  if (sorted) {
    std::sort(resolutions.begin(), resolutions.end());
  }
  return resolutions;
}

namespace internal {
inline std::vector<balancing::Method> avail_normalizations_union(const std::filesystem::path& path,
                                                                 MatrixType matrix_type,
                                                                 MatrixUnit matrix_unit) {
  hic::internal::HiCFileReader reader(path.string());

  const auto& resolutions = reader.header().resolutions;
  if (resolutions.empty()) {
    return {};
  }

  if (resolutions.size() == 1) {
    auto norms = reader.list_avail_normalizations(matrix_type, matrix_unit, resolutions.front());
    std::sort(norms.begin(), norms.end());
    return norms;
  }

  phmap::flat_hash_set<balancing::Method> norms;
  for (const auto& res : resolutions) {
    for (const auto& norm : reader.list_avail_normalizations(matrix_type, matrix_unit, res)) {
      norms.emplace(norm);
    }
  }

  std::vector norms_sorted(norms.begin(), norms.end());
  std::sort(norms_sorted.begin(), norms_sorted.end());
  return norms_sorted;
}

inline std::vector<balancing::Method> avail_normalizations_intersection(
    const std::filesystem::path& path, MatrixType matrix_type, MatrixUnit matrix_unit) {
  hic::internal::HiCFileReader reader(path.string());

  const auto& resolutions = reader.header().resolutions;
  if (resolutions.empty()) {
    return {};
  }

  if (resolutions.size() == 1) {
    auto norms = reader.list_avail_normalizations(matrix_type, matrix_unit, resolutions.front());
    std::sort(norms.begin(), norms.end());
    return norms;
  }

  phmap::flat_hash_map<balancing::Method, std::uint32_t> norms;
  for (const auto& res : resolutions) {
    for (const auto& norm : reader.list_avail_normalizations(matrix_type, matrix_unit, res)) {
      auto [it, inserted] = norms.try_emplace(norm, std::uint32_t{1});
      if (!inserted) {
        it->second++;
      }
    }
  }

  std::vector<balancing::Method> filtered_norms{};
  for (const auto& [norm, count] : norms) {
    if (count == resolutions.size()) {
      filtered_norms.emplace_back(norm);
    }
  }

  std::sort(filtered_norms.begin(), filtered_norms.end());
  return filtered_norms;
}

}  // namespace internal

inline std::vector<balancing::Method> list_normalizations(const std::filesystem::path& path,
                                                          std::string_view policy,
                                                          MatrixType matrix_type,
                                                          MatrixUnit matrix_unit) {
  if (policy == "union") {
    return internal::avail_normalizations_union(path, matrix_type, matrix_unit);
  }
  if (policy == "intersection") {
    return internal::avail_normalizations_intersection(path, matrix_type, matrix_unit);
  }
  throw std::invalid_argument(R"(policy should be either "union" or "intersection")");
}

}  // namespace hictk::hic::utils
