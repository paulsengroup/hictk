// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string_view>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/reference.hpp"

namespace hictk::cooler::utils {

/// Iterable of hictk::File or strings
template <typename N, typename Str>
void merge(Str first_uri, Str last_uri, std::string_view dest_uri, bool overwrite_if_exists = false,
           std::size_t chunk_size = 500'000, std::size_t update_frequency = 10'000'000);

template <typename PixelIt>
void merge(const std::vector<PixelIt>& heads, const std::vector<PixelIt>& tails,
           const BinTable& bins, std::string_view dest_uri, bool overwrite_if_exists = false,
           std::size_t chunk_size = 500'000, std::size_t update_frequency = 10'000'000);

[[nodiscard]] bool equal(std::string_view uri1, std::string_view uri2,
                         bool ignore_attributes = true);
[[nodiscard]] bool equal(const File& clr1, const File& clr2, bool ignore_attributes = true);

[[nodiscard]] std::vector<std::uint32_t> list_resolutions(const std::filesystem::path& path,
                                                          bool sorted = true);

void copy(std::string_view uri1, std::string_view uri2);
void copy(std::string_view uri1, RootGroup dest);

}  // namespace hictk::cooler::utils

#include "./impl/utils_copy_impl.hpp"   // NOLINT
#include "./impl/utils_equal_impl.hpp"  // NOLINT
#include "./impl/utils_impl.hpp"        // NOLINT
#include "./impl/utils_merge_impl.hpp"  // NOLINT
