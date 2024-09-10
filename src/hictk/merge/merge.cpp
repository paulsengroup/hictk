// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <chrono>
#include <cstdint>
#include <filesystem>

#include "hictk/file.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

static void merge_to_cool(const MergeConfig& c) {
  SPDLOG_INFO(FMT_STRING("begin merging {} files into one .{} file..."), c.input_files.size(),
              c.output_format);
  if (c.count_type == "int") {
    utils::merge_to_cool<std::int32_t>(c.input_files.begin(), c.input_files.end(),
                                       c.output_file.string(), c.resolution, c.force, c.chunk_size,
                                       10'000'000, c.compression_lvl);
    return;
  }
  assert(c.count_type == "float");
  utils::merge_to_cool<double>(c.input_files.begin(), c.input_files.end(), c.output_file.string(),
                               c.resolution, c.force, c.chunk_size, 10'000'000, c.compression_lvl);
}

static void merge_to_hic(const MergeConfig& c) {
  SPDLOG_INFO(FMT_STRING("begin merging {} files into one .{} file..."), c.input_files.size(),
              c.output_format);
  utils::merge_to_hic(c.input_files.begin(), c.input_files.end(), c.output_file.string(),
                      c.resolution, c.tmp_dir, c.force, c.chunk_size, c.threads, c.compression_lvl,
                      c.skip_all_vs_all_matrix);
}

int merge_subcmd(const MergeConfig& c) {
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
  SPDLOG_INFO(FMT_STRING("{} size: {:.2f} MB"), c.output_file,
              static_cast<double>(std::filesystem::file_size(c.output_file)) / 1.0e6);

  return 0;
}

}  // namespace hictk::tools
