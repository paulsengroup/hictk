// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <vector>

#include "hictk/hic/file_writer.hpp"
#include "hictk/hic/filestream.hpp"

namespace hictk::hic::internal {

class HiCFileZoomify {
  std::string _path_to_input_hic{};

  HiCFileWriter _hfw{};

 public:
  HiCFileZoomify(std::string_view input_hic, std::string_view output_hic,
                 const std::vector<std::uint32_t>& resolutions, std::size_t n_threads = 1,
                 std::size_t chunk_size = 10'000'000,
                 const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
                 std::uint32_t compression_lvl = 11, bool skip_all_vs_all_matrix = false);
  void zoomify();

 private:
  [[nodiscard]] static HiCFileWriter init_writer(std::string_view input_hic,
                                                 std::string_view output_hic,
                                                 const std::vector<std::uint32_t>& resolution,
                                                 std::size_t n_threads, std::size_t chunk_size,
                                                 const std::filesystem::path& tmpdir,
                                                 std::uint32_t compression_lvl,
                                                 bool skip_all_vs_all_matrix);
  void init();

  [[nodiscard]] std::uint32_t compute_base_resolution(std::uint32_t tgt_resolution) const;

  void ingest_interactions(std::uint32_t resolution);
  void coarsen_interactions(std::uint32_t resolution, std::uint32_t base_resolution);
};

}  // namespace hictk::hic::internal
#include "./impl/file_zoomify_impl.hpp"  // NOLINT
