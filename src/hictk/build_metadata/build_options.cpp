// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/build_options.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <exception>
#include <nlohmann/json.hpp>
#include <string>

namespace hictk::tools {

nlohmann::json get_build_options_json() noexcept {
  try {
    nlohmann::json opts;

    opts.emplace("arch", build_arch());
    opts.emplace("compiler_name", compiler_name());
    opts.emplace("compiler_version", compiler_version());
    opts.emplace("os_name", build_os_name());
    opts.emplace("os_version", build_os_version());
    opts.emplace("build_type", build_type());

    return opts;
  } catch (const std::exception& e) {
    SPDLOG_WARN(FMT_STRING("failed to collect build options: {}"), e.what());
  } catch (...) {
    SPDLOG_WARN("failed to collect build options: unknown error");
  }

  return {};
}

std::string get_build_options(bool pretty) noexcept {
  try {
    if (pretty) {
      return get_build_options_json().dump(2, ' ');
    }
    return get_build_options_json().dump();
  } catch (const std::exception& e) {
    SPDLOG_WARN(FMT_STRING("failed to collect build options: {}"), e.what());
  } catch (...) {
    SPDLOG_WARN("failed to collect build options: unknown error");
  }

  return "{}";
}

}  // namespace hictk::tools
