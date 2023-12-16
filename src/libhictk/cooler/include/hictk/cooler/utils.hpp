// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"

namespace hictk::cooler::utils {

/// Iterable of hictk::File or strings
template <typename N, typename Str>
void merge(Str first_file, Str last_file, std::string_view dest_uri,
           bool overwrite_if_exists = false, std::size_t chunk_size = 500'000,
           std::size_t update_frequency = 10'000'000);

template <typename PixelIt>
void merge(const std::vector<PixelIt>& heads, const std::vector<PixelIt>& tails,
           const BinTable& bins, std::string_view dest_uri, bool overwrite_if_exists = false,
           std::size_t chunk_size = 500'000, std::size_t update_frequency = 10'000'000);

[[nodiscard]] bool equal(std::string_view uri1, std::string_view uri2,
                         bool ignore_attributes = true);
[[nodiscard]] bool equal(const File& clr1, const File& clr2, bool ignore_attributes = true);

[[nodiscard]] std::vector<std::uint32_t> list_resolutions(const std::filesystem::path& path,
                                                          bool sorted = true);

void copy(std::string_view uri1, std::string_view uri2, bool force_overwrite);
void copy(std::string_view uri1, RootGroup dest);

template <typename It>
void rename_chromosomes(std::string_view uri, It first_mapping, It last_mapping);

template <typename NameMap, typename = std::enable_if_t<is_map_v<NameMap>>>
void rename_chromosomes(std::string_view uri, const NameMap& mappings);

template <typename NameMap, typename = std::enable_if_t<is_map_v<NameMap>>>
inline void rename_chromosomes(cooler::Dataset& chrom_dset, const NameMap& mappings);

}  // namespace hictk::cooler::utils

#include "./impl/utils_copy_impl.hpp"
#include "./impl/utils_equal_impl.hpp"
#include "./impl/utils_impl.hpp"
#include "./impl/utils_merge_impl.hpp"
#include "./impl/utils_rename_chroms_impl.hpp"
