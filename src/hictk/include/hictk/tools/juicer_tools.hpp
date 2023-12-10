// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

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

[[nodiscard]] inline std::vector<std::string> generate_juicer_tools_pre_args(
    const ConvertConfig& c, const std::filesystem::path& path_to_pixels,
    const std::filesystem::path& path_to_chrom_sizes, std::size_t processes) {
  assert(processes != 0);
  return {fmt::format(FMT_STRING("-Xmx{}M"), c.juicer_tools_xmx / 1'000'000),
          "-jar",
          c.juicer_tools_jar.string(),
          "pre",
          "-j",
          fmt::to_string(processes),
          "-t",
          c.tmp_dir.string(),
          "-n",
          "-r",
          fmt::format(FMT_STRING("{}"), fmt::join(c.resolutions, ",")),
          path_to_pixels.string(),
          c.path_to_output.string(),
          path_to_chrom_sizes.string()};
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

[[nodiscard]] inline std::unique_ptr<boost::process::child> run_juicer_tools_pre(
    const ConvertConfig& c, const std::filesystem::path& chrom_sizes,
    const std::filesystem::path& pixels, std::size_t processes) {
  const auto cmd = generate_juicer_tools_pre_args(c, pixels, chrom_sizes, processes);
  return std::make_unique<boost::process::child>(find_java().string(), cmd);
}

[[nodiscard]] inline std::unique_ptr<boost::process::child> run_juicer_tools_add_norm(
    const std::filesystem::path& juicer_tools_jar, const std::filesystem::path& path_to_weights,
    const std::filesystem::path& path_to_output, std::size_t juicer_tools_xmx) {
  const auto cmd = generate_juicer_tools_add_norm_args(juicer_tools_jar, path_to_weights,
                                                       path_to_output, juicer_tools_xmx);
  return std::make_unique<boost::process::child>(find_java().string(), cmd);
}

}  // namespace hictk::tools
