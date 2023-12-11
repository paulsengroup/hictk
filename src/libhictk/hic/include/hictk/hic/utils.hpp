// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <cstdint>
#include <filesystem>
#include <vector>

namespace hictk::hic::utils {
[[nodiscard]] bool is_hic_file(const std::filesystem::path &path);
[[nodiscard]] std::vector<std::uint32_t> list_resolutions(const std::filesystem::path &path,
                                                          bool sorted = true);
}  // namespace hictk::hic::utils

#include "./impl/utils_impl.hpp"
