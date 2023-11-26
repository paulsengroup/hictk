// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cassert>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <string_view>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/cooler/cooler.hpp"

namespace hictk::cooler::utils {

namespace internal {
[[nodiscard]] inline std::vector<std::string> get_chrom_names(const cooler::File& clr) {
  std::vector<std::string> names(clr.chromosomes().size());

  std::transform(clr.chromosomes().begin(), clr.chromosomes().end(), names.begin(),
                 [](const hictk::Chromosome& chrom) { return std::string{chrom.name()}; });

  return names;
}

template <typename NameMap>
[[nodiscard]] inline std::vector<std::string>& rename_chromosomes(std::vector<std::string>&& names,
                                                                  const NameMap& mappings) {
  for (auto& name : names) {
    auto it = mappings.find(name);
    if (it != mappings.end()) {
      name = it->second;
    }
  }

  return names;
}

[[nodiscard]] inline std::string find_chrom_with_longest_name(
    const std::vector<std::string>& names) {
  assert(!names.empty());
  return *std::max_element(names.begin(), names.end(), [&](const auto& name1, const auto& name2) {
    return name1.size() < name2.size();
  });
}

}  // namespace internal

template <typename It>
inline void rename_chromosomes(std::string_view uri, It first_mapping, It last_mapping) {
  return rename_chromosomes(
      uri, phmap::btree_map<std::string, std::string>{first_mapping, last_mapping});
}

template <typename NameMap, typename>
inline void rename_chromosomes(std::string_view uri, const NameMap& mappings) {
  if (mappings.empty()) {
    return;
  }
  cooler::File clr(uri);
  auto names = internal::get_chrom_names(clr);
  const auto file_path = clr.path();
  const auto chrom_dset = fmt::format(FMT_STRING("{}/chroms/name"), clr.hdf5_path());
  clr.close();

  names = internal::rename_chromosomes(std::move(names), mappings);

  // NOLINTNEXTLINE(misc-const-correctness)
  HighFive::File h5f(file_path, HighFive::File::ReadWrite);
  const cooler::RootGroup root_grp{h5f.getGroup("/")};
  const auto aprop = h5f.getDataSet(chrom_dset).getAccessPropertyList();

  h5f.unlink(chrom_dset);
  cooler::Dataset dset{root_grp, chrom_dset, internal::find_chrom_with_longest_name(names),
                       HighFive::DataSpace::UNLIMITED, aprop};

  try {
    dset.write(names.begin(), names.end(), 0, true, [&](const auto& name) { return name; });
  } catch (const HighFive::Exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write {} chromosome name(s) to \"{}\": {}"), names.size(),
                    dset.uri(), e.what()));
  }
  assert(dset.size() == names.size());
}
}  // namespace hictk::cooler::utils
