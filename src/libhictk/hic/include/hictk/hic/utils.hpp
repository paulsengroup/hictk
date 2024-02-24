// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <vector>

namespace hictk::hic::utils {

/// Iterable of hictk::hic::File or strings
template <typename Str>
void merge(Str first_file, Str last_file, std::string_view dest_file, std::uint32_t resolution,
           const std::filesystem::path& tmp_dir = std::filesystem::temp_directory_path(),
           bool overwrite_if_exists = false, std::size_t chunk_size = 500'000,
           std::size_t n_threads = 1, std::uint32_t compression_lvl = 11,
           bool skip_all_vs_all = false);

template <typename PixelIt>
void merge(const std::vector<PixelIt>& heads, const std::vector<PixelIt>& tails,
           const BinTable& bins, std::string_view dest_file, std::string_view assembly = "unknown",
           const std::filesystem::path& tmp_dir = std::filesystem::temp_directory_path(),
           bool overwrite_if_exists = false, std::size_t chunk_size = 500'000,
           std::size_t n_threads = 1, std::uint32_t compression_lvl = 11,
           bool skip_all_vs_all = false);

[[nodiscard]] std::vector<std::uint32_t> list_resolutions(const std::filesystem::path& path,
                                                          bool sorted = true);
}  // namespace hictk::hic::utils

#include "./impl/utils_impl.hpp"        // NOLINT
#include "./impl/utils_merge_impl.hpp"  // NOLINT
