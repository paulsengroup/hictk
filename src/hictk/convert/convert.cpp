// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <filesystem>

#include "./common.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

int convert_subcmd(const ConvertConfig& c) {
  auto t0 = std::chrono::steady_clock::now();
  SPDLOG_INFO(FMT_STRING("Converting {} to {} ({} -> {})..."), c.path_to_input, c.path_to_output,
              c.input_format, c.output_format);
  if (c.input_format == "hic") {
    hic_to_cool(c);
  } else {
    cool_to_hic(c);
  }
  auto t1 = std::chrono::steady_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
      1000.0;
  const auto path_to_input = cooler::parse_cooler_uri(c.path_to_input.string()).file_path;
  const auto path_to_output = cooler::parse_cooler_uri(c.path_to_output.string()).file_path;
  SPDLOG_INFO(FMT_STRING("DONE! Processed {} resolution(s) in {:.2f}s!"), c.resolutions.size(),
              delta);
  SPDLOG_INFO(FMT_STRING("{} size: {:.2f} MB"), path_to_input,
              static_cast<double>(std::filesystem::file_size(path_to_input)) / 1.0e6);
  SPDLOG_INFO(FMT_STRING("{} size: {:.2f} MB"), path_to_output,
              static_cast<double>(std::filesystem::file_size(path_to_output)) / 1.0e6);

  return 0;
}
}  // namespace hictk::tools
