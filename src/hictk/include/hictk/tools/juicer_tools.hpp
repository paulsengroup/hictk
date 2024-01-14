// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <boost/process/child.hpp>
#include <boost/process/search_path.hpp>
#include <cassert>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "hictk/tools/config.hpp"

namespace hictk::tools {

[[nodiscard]] inline std::filesystem::path find_java() {
  auto java = boost::process::search_path("java");
  if (java.empty()) {
    throw std::runtime_error("unable to find java in your PATH");
  }
  return java.string();
}

[[nodiscard]] inline std::vector<std::string> generate_juicer_tools_add_norm_args(
    const std::filesystem::path& juicer_tools_jar, const std::filesystem::path& path_to_weights,
    const std::filesystem::path& path_to_output, std::size_t juicer_tools_xmx) {
  return {fmt::format(FMT_STRING("-Xmx{}M"), juicer_tools_xmx / 1'000'000),
          "-jar",
          juicer_tools_jar.string(),
          "addNorm",
          "-j",
          "1",
          path_to_output.string(),
          path_to_weights.string()};
}

[[nodiscard]] inline std::unique_ptr<boost::process::child> run_juicer_tools_add_norm(
    const std::filesystem::path& juicer_tools_jar, const std::filesystem::path& path_to_weights,
    const std::filesystem::path& path_to_output, std::size_t juicer_tools_xmx) {
  const auto cmd = generate_juicer_tools_add_norm_args(juicer_tools_jar, path_to_weights,
                                                       path_to_output, juicer_tools_xmx);
  return std::make_unique<boost::process::child>(find_java().string(), cmd);
}

}  // namespace hictk::tools
