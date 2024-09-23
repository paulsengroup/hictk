// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cstdint>
#include <string_view>
#include <vector>

#include "./zoomify.hpp"
#include "hictk/hic/file_zoomify.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

static void print_zooming_plan_hic(std::string_view path_to_input,
                                   const std::vector<std::uint32_t>& resolutions) {
  const auto avail_resolutions = hic::utils::list_resolutions(path_to_input);
  for (const auto& res : resolutions) {
    const auto match = std::find(avail_resolutions.begin(), avail_resolutions.end(), res);
    if (match != avail_resolutions.end()) {
      SPDLOG_INFO(FMT_STRING("copying resolution {} from \"{}\""), res, path_to_input);
    } else {
      auto base_resolution = avail_resolutions.front();
      for (const auto& avail_res : resolutions) {
        if (avail_res >= res) {
          break;
        }
        if (res % avail_res == 0) {
          base_resolution = avail_res;
        }
      }
      SPDLOG_INFO(FMT_STRING("generating {} resolution from {} ({}x)"), res, base_resolution,
                  res / base_resolution);
    }
  }
}

void zoomify_hic(const ZoomifyConfig& c) {
  if (c.force) {
    std::filesystem::remove(c.path_to_output);
  }

  print_zooming_plan_hic(c.path_to_input.string(), c.resolutions);

  const internal::TmpDir tmpdir{c.tmp_dir, true};
  hic::internal::HiCFileZoomify{c.path_to_input.string(),
                                c.path_to_output.string(),
                                c.resolutions,
                                c.threads,
                                c.batch_size,
                                tmpdir(),
                                c.compression_lvl,
                                c.skip_all_vs_all_matrix}
      .zoomify();
}

}  // namespace hictk::tools
