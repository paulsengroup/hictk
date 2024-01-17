// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <chrono>
#include <cstdint>
#include <filesystem>

#include "hictk/cooler/utils.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

static void merge_coolers(const MergeConfig& c) {
  SPDLOG_INFO(FMT_STRING("begin merging {} coolers..."), c.input_files.size());
  cooler::utils::merge<std::int32_t>(c.input_files.begin(), c.input_files.end(),
                                     c.output_file.string(), c.force, c.chunk_size, 10'000'000,
                                     c.compression_lvl);
}

static void merge_hics(const MergeConfig& c) {
  SPDLOG_INFO(FMT_STRING("begin merging {} .hic files..."), c.input_files.size());
  const internal::TmpDir tmpdir{c.tmp_dir};
  hic::utils::merge(c.input_files.begin(), c.input_files.end(), c.output_file.string(),
                    c.resolution, tmpdir(), c.force, c.chunk_size, c.threads, c.compression_lvl);
}

int merge_subcmd(const MergeConfig& c) {
  const auto t0 = std::chrono::system_clock::now();
  if (c.output_format == "cool") {
    merge_coolers(c);
  } else {
    merge_hics(c);
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
