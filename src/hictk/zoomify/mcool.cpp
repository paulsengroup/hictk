// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <variant>
#include <vector>

#include "./zoomify.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"

namespace hictk::tools {

void zoomify_many_cooler(std::string_view in_uri, std::string_view out_path,
                         const std::vector<std::uint32_t>& resolutions, bool copy_base_resolution,
                         bool force, std::uint32_t compression_lvl) {
  const cooler::File clr(in_uri);
  auto mclr = cooler::MultiResFile::create(out_path, cooler::File(in_uri).chromosomes(), force);

  SPDLOG_INFO(FMT_STRING("coarsening cooler at {} {} times ({} -> {})"), clr.uri(),
              resolutions.size(), clr.resolution(), fmt::join(resolutions, " -> "));

  if (copy_base_resolution) {
    assert(resolutions.front() == clr.resolution());
    mclr.copy_resolution(clr);
  } else {
    assert(resolutions.size() > 1);
    zoomify_once_cooler(cooler::File(in_uri), mclr.init_resolution(resolutions[1]), resolutions[1],
                        compression_lvl);
  }

  std::visit(
      [&]([[maybe_unused]] auto count_type) {
        using PixelT = decltype(count_type);
        for (std::size_t i = 1; i < resolutions.size(); ++i) {
          mclr.create_resolution<PixelT>(resolutions[i]);
        }
      },
      clr.pixel_variant());
}

}  // namespace hictk::tools
