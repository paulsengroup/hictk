// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "./common.hpp"
#include "./load_pairs.hpp"
#include "./load_pixels.hpp"
#include "hictk/hic/file_writer.hpp"

namespace hictk::tools {

static Stats ingest_pixels_hic(std::string_view uri, const std::filesystem::path& tmp_dir,
                               const Reference& chromosomes, std::uint32_t bin_size,
                               const std::string& assembly, std::int64_t offset,
                               bool skip_all_vs_all_matrix, Format format, std::size_t threads,
                               std::size_t batch_size, std::uint32_t compression_lvl, bool force) {
  SPDLOG_INFO(FMT_STRING("begin loading pixels into a .hic file..."));

  if (force) {
    std::filesystem::remove(uri);
  }

  hic::internal::HiCFileWriter hf(uri, chromosomes, {bin_size}, assembly, threads, batch_size,
                                  tmp_dir, compression_lvl, skip_all_vs_all_matrix);

  std::vector<ThinPixel<float>> write_buffer(batch_size);
  return ingest_pixels(std::move(hf), write_buffer, format, offset);
}

inline Stats ingest_pairs_hic(std::string_view uri, const std::filesystem::path& tmp_dir,
                              const Reference& chromosomes, std::uint32_t bin_size,
                              const std::string& assembly, std::int64_t offset,
                              bool skip_all_vs_all_matrix, Format format, std::size_t threads,
                              std::size_t batch_size, std::uint32_t compression_lvl, bool force) {
  if (force) {
    std::filesystem::remove(uri);
  }

  hic::internal::HiCFileWriter hf(uri, chromosomes, {bin_size}, assembly, threads, batch_size,
                                  tmp_dir, compression_lvl, skip_all_vs_all_matrix);

  std::vector<ThinPixel<float>> buffer(batch_size);
  return ingest_pairs(std::move(hf), buffer, format, offset);
}

}  // namespace hictk::tools
