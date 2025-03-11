// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string_view>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk::hic::utils {

// NOLINTBEGIN(*-avoid-magic-numbers)

/// Iterable of hictk::hic::File or strings
template <typename Str>
void merge(
    Str first_file, Str last_file, std::string_view dest_file, std::uint32_t resolution,
    const std::filesystem::path& tmp_dir = hictk::internal::TmpDir::default_temp_directory_path(),
    bool overwrite_if_exists = false, std::size_t chunk_size = 500'000, std::size_t n_threads = 1,
    std::uint32_t compression_lvl = 11, bool skip_all_vs_all = false);

template <typename PixelIt>
void merge(
    const std::vector<PixelIt>& heads, const std::vector<PixelIt>& tails, const BinTable& bins,
    std::string_view dest_file, std::string_view assembly = "unknown",
    const std::filesystem::path& tmp_dir = hictk::internal::TmpDir::default_temp_directory_path(),
    bool overwrite_if_exists = false, std::size_t chunk_size = 500'000, std::size_t n_threads = 1,
    std::uint32_t compression_lvl = 11, bool skip_all_vs_all = false);

// NOLINTEND(*-avoid-magic-numbers)

[[nodiscard]] std::vector<std::uint32_t> list_resolutions(const std::filesystem::path& path,
                                                          bool sorted = true);
[[nodiscard]] std::vector<balancing::Method> list_normalizations(
    const std::filesystem::path& path, std::string_view policy = "union",
    MatrixType matrix_type = MatrixType::observed, MatrixUnit matrix_unit = MatrixUnit::BP);

}  // namespace hictk::hic::utils

#include "./impl/utils_impl.hpp"        // NOLINT
#include "./impl/utils_merge_impl.hpp"  // NOLINT
