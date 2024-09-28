// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

void zoomify_hic(const ZoomifyConfig& c);

void zoomify_cooler(const ZoomifyConfig& c, bool output_is_multires);

void zoomify_once_cooler(const cooler::File& clr1, cooler::RootGroup entrypoint2,
                         std::uint32_t resolution, std::uint32_t compression_lvl);

void zoomify_many_cooler(std::string_view in_uri, std::string_view out_path,
                         const std::vector<std::uint32_t>& resolutions, bool copy_base_resolution,
                         bool force, std::uint32_t compression_lvl);

}  // namespace hictk::tools
