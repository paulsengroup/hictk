// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./merge.hpp"

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <filesystem>

#include "hictk/tools/config.hpp"

namespace hictk::tools {

int merge_subcmd(const MergeConfig& c) {  // NOLINT(misc-use-internal-linkage)
  const auto t0 = std::chrono::system_clock::now();
  if (c.output_format == "cool") {
    merge_to_cool(c);
  } else {
    assert(c.output_format == "hic");
    merge_to_hic(c);
  }

  const auto t1 = std::chrono::system_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
      1000.0;

  SPDLOG_INFO(FMT_STRING("DONE! Merging {} files took {:.2f}s!"), c.input_files.size(), delta);
  // NOLINTBEGIN(*-avoid-magic-numbers)
  SPDLOG_INFO(FMT_STRING("{} size: {:.2f} MB"), c.output_file,
              static_cast<double>(std::filesystem::file_size(c.output_file)) / 1.0e6);
  // NOLINTEND(*-avoid-magic-numbers)

  return 0;
}

}  // namespace hictk::tools
