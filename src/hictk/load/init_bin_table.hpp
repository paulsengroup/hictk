// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>

#include "hictk/bin_table.hpp"

namespace hictk::tools {

[[nodiscard]] BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes,
                                      std::uint32_t bin_size);

[[nodiscard]] BinTable init_bin_table(const std::filesystem::path& path_to_bin_table);

[[nodiscard]] BinTable init_bin_table(const std::filesystem::path& path_to_chrom_sizes,
                                      const std::filesystem::path& path_to_bin_table,
                                      std::uint32_t bin_size);

}  // namespace hictk::tools
