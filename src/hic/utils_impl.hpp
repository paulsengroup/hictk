// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <string>

namespace hictk::hic::utils {
inline bool is_hic_file(const std::filesystem::path& path) {
  return internal::HiCFileReader::checkMagicString(path.string());
}

inline std::vector<std::uint32_t> list_resolutions(const std::filesystem::path& path) {
  return hic::internal::HiCFileReader(path.string()).header().resolutions;
}
}  // namespace hictk::hic::utils
