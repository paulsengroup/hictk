// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>

#include "./common.hpp"
#include "hictk/reference.hpp"

namespace hictk::tools {

[[nodiscard]] Stats ingest_pixels_hic(PixelQueue<float>& pixel_queue,
                                      const std::atomic<bool>& early_return, std::string_view uri,
                                      const std::filesystem::path& tmp_dir,
                                      const Reference& chromosomes, std::uint32_t bin_size,
                                      const std::string& assembly, bool skip_all_vs_all_matrix,
                                      std::size_t threads, std::size_t batch_size,
                                      std::uint32_t compression_lvl, bool validate, bool force);

[[nodiscard]] Stats ingest_pairs_hic(PixelQueue<float>& pixel_queue,
                                     const std::atomic<bool>& early_return, std::string_view uri,
                                     const std::filesystem::path& tmp_dir,
                                     const Reference& chromosomes, std::uint32_t bin_size,
                                     const std::string& assembly, bool skip_all_vs_all_matrix,
                                     std::size_t threads, std::size_t batch_size,
                                     std::uint32_t compression_lvl, bool validate, bool force);

}  // namespace hictk::tools
