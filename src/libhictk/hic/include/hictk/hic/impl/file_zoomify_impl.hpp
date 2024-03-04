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
                                      std::uint32_t compression_lvl, bool skip_all_vs_all_matrix)
    : _path_to_input_hic(std::string{input_hic}),
      _hfw(init_writer(input_hic, output_hic, resolutions, n_threads, chunk_size, tmpdir,
                       compression_lvl, skip_all_vs_all_matrix)) {
  init();
}

inline void HiCFileZoomify::init() {
  const auto avail_resolutions = hic::utils::list_resolutions(_path_to_input_hic);
  const auto base_resolution = avail_resolutions.front();
  if (base_resolution != _hfw.resolutions().front()) {
    const auto tgt_resolution = _hfw.resolutions().front();
    coarsen_interactions(tgt_resolution, base_resolution);
  }

  for (const auto& res : _hfw.resolutions()) {
    const auto res_avail = std::find(avail_resolutions.begin(), avail_resolutions.end(), res) !=
                           avail_resolutions.end();
    if (res_avail) {
      ingest_interactions(res);
    }
  }
}

inline void HiCFileZoomify::ingest_interactions(std::uint32_t resolution) {
  // TODO: check if .hic is version 9
  // if it is, copy blocks directly
  // if it isn't copy pixels

  SPDLOG_INFO(FMT_STRING("[{} bp] ingesting interactions..."), resolution);
  const File hf(_path_to_input_hic, resolution);
  const auto sel = hf.fetch();
  _hfw.add_pixels(resolution, sel.begin<float>(), sel.end<float>());
}

inline void HiCFileZoomify::coarsen_interactions(std::uint32_t resolution,
                                                 std::uint32_t base_resolution) {
  if (resolution % base_resolution != 0) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("unable to generate pixels at resolution {} by coarsening resolution {}"),
        resolution, base_resolution));
  }

  SPDLOG_INFO(FMT_STRING("[{} bp] coarsening interactions from res {} ({}x)..."), resolution,
              base_resolution, resolution / base_resolution);
  const File hf(_path_to_input_hic, base_resolution);
  const auto sel1 = hf.fetch();
  const auto sel2 = transformers::CoarsenPixels(sel1.begin<float>(), sel1.end<float>(),
                                                hf.bins_ptr(), resolution / base_resolution);
  _hfw.add_pixels(resolution, sel2.begin(), sel2.end());
}

inline void HiCFileZoomify::zoomify() { _hfw.serialize(); }

inline HiCFileWriter HiCFileZoomify::init_writer(std::string_view input_hic,
                                                 std::string_view output_hic,
                                                 const std::vector<std::uint32_t>& resolutions,
                                                 std::size_t n_threads, std::size_t chunk_size,
                                                 const std::filesystem::path& tmpdir,
                                                 std::uint32_t compression_lvl,
                                                 bool skip_all_vs_all_matrix) {
  auto resolutions_ = resolutions;
  std::sort(resolutions_.begin(), resolutions_.end());

  const auto avail_resolutions = hic::utils::list_resolutions(input_hic);
  const File hf(std::string{input_hic}, avail_resolutions.back());

  return HiCFileWriter{output_hic,    hf.chromosomes(), resolutions_,
                       hf.assembly(), n_threads,        chunk_size,
                       tmpdir,        compression_lvl,  skip_all_vs_all_matrix};
}

}  // namespace hictk::hic::internal
