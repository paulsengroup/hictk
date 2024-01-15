// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <highfive/H5File.hpp>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/pixel.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/transformers/coarsen.hpp"

namespace hictk::tools {

void zoomify_once(const cooler::File& clr1, cooler::RootGroup entrypoint2, std::uint32_t resolution,
                  std::uint32_t compression_lvl) {
  auto attrs = cooler::Attributes::init(clr1.bin_size());
  attrs.assembly = clr1.attributes().assembly;

  auto clr2 = cooler::File::create(std::move(entrypoint2), clr1.chromosomes(), resolution, attrs,
                                   cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl);

  std::vector<ThinPixel<std::int32_t>> buffer{500'000};
  cooler::MultiResFile::coarsen(clr1, clr2, buffer);
}

void zoomify_once(std::string_view uri1, std::string_view uri2, std::uint32_t resolution,
                  bool force, std::uint32_t compression_lvl) {
  const cooler::File clr1(uri1);

  SPDLOG_INFO(FMT_STRING("coarsening cooler at {} once ({} -> {})"), clr1.uri(), clr1.bin_size(),
              resolution);

  auto mode = force ? HighFive::File::Overwrite : HighFive::File::Create;
  cooler::RootGroup entrypoint2{HighFive::File(std::string{uri2}, mode).getGroup("/")};

  return zoomify_once(clr1, std::move(entrypoint2), resolution, compression_lvl);
}

void zoomify_many(std::string_view in_uri, std::string_view out_path,
                  const std::vector<std::uint32_t>& resolutions, bool copy_base_resolution,
                  bool force, std::uint32_t compression_lvl) {
  const cooler::File clr(in_uri);
  auto mclr = cooler::MultiResFile::create(out_path, cooler::File(in_uri).chromosomes(), force);

  SPDLOG_INFO(FMT_STRING("coarsening cooler at {} {} times ({} -> {})"), clr.uri(),
              resolutions.size(), clr.bin_size(), fmt::join(resolutions, " -> "));

  if (copy_base_resolution) {
    assert(resolutions.front() == clr.bin_size());
    mclr.copy_resolution(clr);
  } else {
    assert(resolutions.size() > 1);
    zoomify_once(cooler::File(in_uri), mclr.init_resolution(resolutions[1]), resolutions[1],
                 compression_lvl);
  }

  for (std::size_t i = 1; i < resolutions.size(); ++i) {
    mclr.create_resolution(resolutions[i]);
  }
}

int zoomify_subcmd(const ZoomifyConfig& c) {
  const auto t0 = std::chrono::system_clock::now();
  const auto output_is_multires = c.copy_base_resolution || c.resolutions.size() > 2;

  if (output_is_multires) {
    zoomify_many(c.input_uri, c.output_path, c.resolutions, c.copy_base_resolution, c.force,
                 c.compression_lvl);
  } else {
    zoomify_once(c.input_uri, c.output_path, c.resolutions.back(), c.force, c.compression_lvl);
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
