// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./zoomify.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/hic/file_zoomify.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/zoomify.hpp"

namespace hictk::tools {

int run_subcmd(const ZoomifyConfig& c) {  // NOLINT(misc-use-internal-linkage)
  const auto output_is_multires = c.copy_base_resolution || c.resolutions.size() > 2;
  const auto t0 = std::chrono::system_clock::now();
  if (c.output_format == "hic") {
    zoomify_hic(c);
  } else {
    zoomify_cooler(c, output_is_multires);
  }
  const auto t1 = std::chrono::system_clock::now();
  const auto delta =
      static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
      1000.0;
  SPDLOG_INFO(FMT_STRING("DONE! Processed {} resolution(s) in {:.2f}s!"),
              c.resolutions.size() - !output_is_multires, delta);

  return 0;
}
}  // namespace hictk::tools
