// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "./merge.hpp"
#include "hictk/file.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

// NOLINTNEXTLINE(misc-use-internal-linkage)
void merge_to_hic(const MergeConfig& c) {
  SPDLOG_INFO(FMT_STRING("begin merging {} files into one .{} file..."), c.input_files.size(),
              c.output_format);
  utils::merge_to_hic(c.input_files.begin(), c.input_files.end(), c.output_file.string(),
                      c.resolution, c.tmp_dir, c.force, c.chunk_size, c.threads, c.compression_lvl,
                      c.skip_all_vs_all_matrix);
}

}  // namespace hictk::tools
