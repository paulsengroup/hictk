// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstddef>
#include <cstdint>

#include "hictk/file.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

// NOLINTNEXTLINE(misc-use-internal-linkage)
void merge_to_cool(const MergeConfig& c) {
  constexpr std::size_t update_freq = 10'000'000;
  SPDLOG_INFO(FMT_STRING("begin merging {} files into one .{} file..."), c.input_files.size(),
              c.output_format);
  if (c.count_type == "int") {
    utils::merge_to_cool<std::int32_t>(c.input_files.begin(), c.input_files.end(),
                                       c.output_file.string(), c.resolution, c.force, c.chunk_size,
                                       update_freq, c.compression_lvl);
    return;
  }
  assert(c.count_type == "float");
  utils::merge_to_cool<double>(c.input_files.begin(), c.input_files.end(), c.output_file.string(),
                               c.resolution, c.force, c.chunk_size, update_freq, c.compression_lvl);
}

}  // namespace hictk::tools
