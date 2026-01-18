// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/validation.hpp"

#include <filesystem>

#include "hictk/hic/file_reader.hpp"

namespace hictk::hic::utils {
bool is_hic_file(const std::filesystem::path& path) {
  return internal::HiCFileReader::checkMagicString(path.string());
}

}  // namespace hictk::hic::utils
