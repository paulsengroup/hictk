// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>

namespace hictk::hic::utils {
[[nodiscard]] bool is_hic_file(const std::filesystem::path &path);
}  // namespace hictk::hic::utils

#include "./impl/validation_impl.hpp"  // NOLINT
