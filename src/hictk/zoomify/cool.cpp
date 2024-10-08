// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cstdint>
#include <highfive/H5File.hpp>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

#include "./zoomify.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

void zoomify_once_cooler(const cooler::File& clr1, cooler::RootGroup entrypoint2,
                         std::uint32_t resolution, std::uint32_t compression_lvl) {
  auto attrs = cooler::Attributes::init(clr1.resolution());
  attrs.assembly = clr1.attributes().assembly;

  std::visit(
      [&]([[maybe_unused]] auto count_type) {
        using PixelT = decltype(count_type);
        auto clr2 = cooler::File::create<PixelT>(
            std::move(entrypoint2), clr1.chromosomes(), resolution, attrs,
            cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl);

        std::vector<ThinPixel<PixelT>> buffer{500'000};
        cooler::MultiResFile::coarsen(clr1, clr2, buffer);
      },
      clr1.pixel_variant());
}

static void zoomify_once_cooler(std::string_view uri1, std::string_view uri2,
                                std::uint32_t resolution, bool force,
                                std::uint32_t compression_lvl) {
  const cooler::File clr1(uri1);

  SPDLOG_INFO(FMT_STRING("coarsening cooler at {} once ({} -> {})"), clr1.uri(), clr1.resolution(),
              resolution);

  auto mode = force ? HighFive::File::Overwrite : HighFive::File::Create;
  cooler::RootGroup entrypoint2{HighFive::File(std::string{uri2}, mode).getGroup("/")};

  zoomify_once_cooler(clr1, std::move(entrypoint2), resolution, compression_lvl);
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

void zoomify_cooler(const ZoomifyConfig& c, bool output_is_multires) {
  if (output_is_multires) {
    zoomify_many_cooler(c.path_to_input.string(), c.path_to_output.string(), c.resolutions,
                        c.copy_base_resolution, c.force, c.compression_lvl);
    return;
  }
  zoomify_once_cooler(c.path_to_input.string(), c.path_to_output.string(), c.resolutions.back(),
                      c.force, c.compression_lvl);
}

}  // namespace hictk::tools
