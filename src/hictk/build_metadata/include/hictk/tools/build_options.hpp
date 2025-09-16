// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <nlohmann/json.hpp>
#include <string>
#include <string_view>

namespace hictk::tools {

[[nodiscard]] nlohmann::json get_build_options_json() noexcept;
[[nodiscard]] std::string get_build_options(bool pretty = false) noexcept;

[[nodiscard]] constexpr std::string_view build_os_name() noexcept {
#ifdef HICTK_SYSTEM_NAME
  return HICTK_SYSTEM_NAME;
#else
  return "unknown";
#endif
}

[[nodiscard]] constexpr std::string_view build_os_version() noexcept {
#ifdef HICTK_SYSTEM_VERSION
  return HICTK_SYSTEM_VERSION;
#else
  return "unknown";
#endif
}

[[nodiscard]] constexpr std::string_view build_arch() noexcept {
#ifdef HICTK_SYSTEM_PROCESSOR
  return HICTK_SYSTEM_PROCESSOR;
#else
  return "unknown";
#endif
}

[[nodiscard]] constexpr std::string_view build_type() noexcept {
#ifdef HICTK_BUILD_TYPE
  return HICTK_BUILD_TYPE;
#else
  return "unknown";
#endif
}

[[nodiscard]] constexpr std::string_view compiler_name() noexcept {
#ifdef HICTK_CXX_COMPILER_ID
  return HICTK_CXX_COMPILER_ID;
#else
  return "unknown";
#endif
}

[[nodiscard]] constexpr std::string_view compiler_version() noexcept {
#ifdef HICTK_CXX_COMPILER_VERSION
  return HICTK_CXX_COMPILER_VERSION;
#else
  return "unknown";
#endif
}

}  // namespace hictk::tools
