// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <string>

namespace hictk::utils {
inline bool is_hic_file(const std::filesystem::path& path) {
  return internal::HiCFileStream::checkMagicString(path.string());
}
}  // namespace hictk::utils
