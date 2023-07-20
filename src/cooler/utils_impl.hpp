// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/std.h>

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <string>
#include <vector>

#include "hictk/cooler.hpp"

namespace hictk::cooler::utils {

inline std::vector<std::uint32_t> list_resolutions(const std::filesystem::path &path, bool sorted) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  try {
    if (!is_multires_file(path.string(), false)) {
      throw std::runtime_error("not a valid .mcool file");
    }

    const HighFive::File fp(path.string(), HighFive::File::ReadOnly);
    auto root_grp = fp.getGroup("/resolutions");

    const auto resolutions_ = root_grp.listObjectNames();
    std::vector<std::uint32_t> resolutions(resolutions_.size());
    std::transform(resolutions_.begin(), resolutions_.end(), resolutions.begin(),
                   [](const auto &res) {
                     return hictk::internal::parse_numeric_or_throw<std::uint32_t>(res);
                   });

    if (sorted) {
      std::sort(resolutions.begin(), resolutions.end());
    }
    return resolutions;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to read resolutions from {}: {}"), path, e.what()));
  }
}

}  // namespace hictk::cooler::utils
