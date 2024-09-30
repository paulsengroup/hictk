// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string>
#include <string_view>
#include <utility>

#include "hictk/tools/toml.hpp"

namespace hictk::tools {

[[nodiscard]] std::pair<int, toml::table> validate_hic(const std::string& path, bool exhaustive);

[[nodiscard]] std::pair<int, toml::table> validate_cooler(std::string_view path,
                                                          bool validate_index);

[[nodiscard]] std::pair<int, toml::table> validate_mcool(std::string_view path, bool validate_index,
                                                         bool exhaustive);

[[nodiscard]] std::pair<int, toml::table> validate_scool(std::string_view path, bool validate_index,
                                                         bool exhaustive);
}  // namespace hictk::tools
