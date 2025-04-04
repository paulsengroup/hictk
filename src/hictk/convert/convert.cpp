// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/convert.hpp"

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <chrono>
#include <filesystem>

#include "./common.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/string_utils.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

int run_subcmd(const ConvertConfig& c) {  // NOLINT(misc-use-internal-linkage)
  auto t0 = std::chrono::steady_clock::now();
  SPDLOG_INFO(FMT_STRING("Converting {} to {} ({} -> {})..."), c.path_to_input, c.path_to_output,
              c.input_format, c.output_format);
  if (c.input_format == "hic") {
    assert(internal::ends_with(c.output_format, "cool"));
    hic_to_cool(c);
  } else {
    assert(internal::starts_with(c.output_format, "hic"));
    cool_to_hic(c);
  }
  auto t1 = std::chrono::steady_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
      1000.0;
  const auto path_to_input = cooler::parse_cooler_uri(c.path_to_input.string()).file_path;
  const auto path_to_output = cooler::parse_cooler_uri(c.path_to_output.string()).file_path;
  // NOLINTBEGIN(*-avoid-magic-numbers)
  SPDLOG_INFO(FMT_STRING("DONE! Processed {} resolution(s) in {:.2f}s!"), c.resolutions.size(),
              delta);
  SPDLOG_INFO(FMT_STRING("{} size: {:.2f} MB"), path_to_input,
              static_cast<double>(std::filesystem::file_size(path_to_input)) / 1.0e6);
  SPDLOG_INFO(FMT_STRING("{} size: {:.2f} MB"), path_to_output,
              static_cast<double>(std::filesystem::file_size(path_to_output)) / 1.0e6);
  // NOLINTEND(*-avoid-magic-numbers)

  return 0;
}
}  // namespace hictk::tools
