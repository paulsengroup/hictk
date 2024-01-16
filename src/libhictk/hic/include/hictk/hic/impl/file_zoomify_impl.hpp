// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <vector>

#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/version.hpp"

namespace hictk::hic::internal {

inline HiCFileZoomify::HiCFileZoomify(std::string_view input_hic, std::string_view output_hic,
                                      const std::vector<std::uint32_t>& resolutions,
                                      std::size_t n_threads, std::size_t chunk_size,
                                      const std::filesystem::path& tmpdir,
                                      std::uint32_t compression_lvl)
    : _path_to_input_hic(std::string{input_hic}),
      _hfw(init_writer(input_hic, output_hic, resolutions, n_threads, chunk_size, tmpdir,
                       compression_lvl)) {
  const auto avail_resolutions = hic::utils::list_resolutions(input_hic);
  const auto base_resolution = avail_resolutions.front();
  for (const auto& res : resolutions) {
    if (res % base_resolution != 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("unable to generate resolution {} from base resolution {}"), res,
                      base_resolution));
    }
  }
  init();
}

inline void HiCFileZoomify::init() {
  const auto avail_resolutions = hic::utils::list_resolutions(_path_to_input_hic);

  // TODO: check if .hic is version 9
  // if it is, copy blocks directly
  // if it isn't copy pixels
  for (const auto& res : _hfw.resolutions()) {
    const auto res_avail = std::find(avail_resolutions.begin(), avail_resolutions.end(), res) !=
                           avail_resolutions.end();
    if (res_avail) {
      const File hf(_path_to_input_hic, res);
      const auto sel = hf.fetch();
      _hfw.add_pixels(res, sel.begin<float>(), sel.end<float>());
    }
  }
}

inline void HiCFileZoomify::zoomify() { _hfw.serialize(); }

inline HiCFileWriter HiCFileZoomify::init_writer(std::string_view input_hic,
                                                 std::string_view output_hic,
                                                 const std::vector<std::uint32_t>& resolutions,
                                                 std::size_t n_threads, std::size_t chunk_size,
                                                 const std::filesystem::path& tmpdir,
                                                 std::uint32_t compression_lvl) {
  auto resolutions_ = resolutions;
  std::sort(resolutions_.begin(), resolutions_.end());

  const auto avail_resolutions = hic::utils::list_resolutions(input_hic);
  const File hf(std::string{input_hic}, avail_resolutions.back());

  return HiCFileWriter{output_hic, hf.chromosomes(), resolutions_, hf.assembly(),
                       n_threads,  chunk_size,       tmpdir,       compression_lvl};
}

}  // namespace hictk::hic::internal
