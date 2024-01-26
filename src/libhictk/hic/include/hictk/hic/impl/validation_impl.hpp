// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "hictk/hic/file_reader.hpp"

namespace hictk::hic::utils {
inline bool is_hic_file(const std::filesystem::path& path) {
  return hictk::hic::internal::HiCFileReader::checkMagicString(path.string());
}

}  // namespace hictk::hic::utils
