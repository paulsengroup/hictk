// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
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
#include "hictk/hic/file_zoomify.hpp"
#include "hictk/pixel.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void zoomify_once_cooler(const cooler::File& clr1, cooler::RootGroup entrypoint2,
                         std::uint32_t resolution, std::uint32_t compression_lvl) {
  auto attrs = cooler::Attributes::init(clr1.bin_size());
  attrs.assembly = clr1.attributes().assembly;

  auto clr2 = cooler::File::create(std::move(entrypoint2), clr1.chromosomes(), resolution, attrs,
                                   cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl);

  std::vector<ThinPixel<std::int32_t>> buffer{500'000};
  cooler::MultiResFile::coarsen(clr1, clr2, buffer);
}

void zoomify_once_cooler(std::string_view uri1, std::string_view uri2, std::uint32_t resolution,
                         bool force, std::uint32_t compression_lvl) {
  const cooler::File clr1(uri1);

  SPDLOG_INFO(FMT_STRING("coarsening cooler at {} once ({} -> {})"), clr1.uri(), clr1.bin_size(),
              resolution);

  auto mode = force ? HighFive::File::Overwrite : HighFive::File::Create;
  cooler::RootGroup entrypoint2{HighFive::File(std::string{uri2}, mode).getGroup("/")};

  return zoomify_once_cooler(clr1, std::move(entrypoint2), resolution, compression_lvl);
}

void zoomify_many_cooler(std::string_view in_uri, std::string_view out_path,
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
    zoomify_once_cooler(cooler::File(in_uri), mclr.init_resolution(resolutions[1]), resolutions[1],
                        compression_lvl);
  }

  for (std::size_t i = 1; i < resolutions.size(); ++i) {
    mclr.create_resolution(resolutions[i]);
  }
}

void print_zooming_plan_hic(std::string_view path_to_input,
                            const std::vector<std::uint32_t>& resolutions) {
  const auto avail_resolutions = hic::utils::list_resolutions(path_to_input);
  for (const auto& res : resolutions) {
    const auto match = std::find(avail_resolutions.begin(), avail_resolutions.end(), res);
    if (match != avail_resolutions.end()) {
      SPDLOG_INFO(FMT_STRING("copying resolution {} from \"{}\""), res, path_to_input);
    } else {
      auto base_resolution = resolutions.front();
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

  const internal::TmpDir tmpdir{c.tmp_dir};
  hic::internal::HiCFileZoomify{c.path_to_input.string(),
                                c.path_to_output.string(),
                                c.resolutions,
                                c.threads,
                                c.batch_size,
                                tmpdir(),
                                c.compression_lvl}
      .zoomify();
}

void zoomify_cooler(const ZoomifyConfig& c, bool output_is_multires) {
  if (output_is_multires) {
    zoomify_many_cooler(c.path_to_input.string(), c.path_to_output.string(), c.resolutions,
                        c.copy_base_resolution, c.force, c.compression_lvl);
    return;
  }
  zoomify_once_cooler(c.path_to_input.string(), c.path_to_output.string(), c.resolutions.back(),
                      c.force, c.compression_lvl);
}

int zoomify_subcmd(const ZoomifyConfig& c) {
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
