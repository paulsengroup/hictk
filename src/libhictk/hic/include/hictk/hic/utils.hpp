// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <vector>

namespace hictk::hic::utils {
[[nodiscard]] std::vector<std::uint32_t> list_resolutions(const std::filesystem::path &path,
                                                          bool sorted = true);
}  // namespace hictk::hic::utils

#include "./impl/utils_impl.hpp"  // NOLINT
