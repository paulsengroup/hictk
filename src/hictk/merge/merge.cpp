// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <cstdint>
#include <filesystem>

#include "hictk/cooler/utils.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

int merge_subcmd(const MergeConfig& c) {
  SPDLOG_INFO(FMT_STRING("begin merging {} coolers..."), c.input_uris.size());
  const auto t0 = std::chrono::system_clock::now();
  cooler::utils::merge<std::int32_t>(c.input_uris.begin(), c.input_uris.end(),
                                     c.output_uri.string(), c.force, c.chunk_size, 10'000'000);
  const auto t1 = std::chrono::system_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
      1000.0;
  SPDLOG_INFO(FMT_STRING("DONE! Merging {} coolers took {:.2f}s!"), c.input_uris.size(), delta);
  SPDLOG_INFO(FMT_STRING("{} size: {:.2f} MB"), c.output_uri,
              static_cast<double>(std::filesystem::file_size(c.output_uri)) / 1.0e6);

  return 0;
}

}  // namespace hictk::tools
