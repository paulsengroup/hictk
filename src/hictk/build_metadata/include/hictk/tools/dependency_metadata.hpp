// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nlohmann/json.hpp>
#include <string>

namespace hictk::tools {

[[nodiscard]] nlohmann::json get_dependency_versions_json() noexcept;
[[nodiscard]] std::string get_dependency_versions(bool pretty = false) noexcept;

}  // namespace hictk::tools
