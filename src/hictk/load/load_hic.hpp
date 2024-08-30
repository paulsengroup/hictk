// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "./common.hpp"
#include "./load_pairs.hpp"
#include "./load_pixels.hpp"
#include "hictk/hic/file_writer.hpp"

namespace hictk::tools {

static Stats ingest_pixels_hic(PixelQueue<float>& pixel_queue,
                               const std::atomic<bool>& early_return, std::string_view uri,
                               const std::filesystem::path& tmp_dir, const Reference& chromosomes,
                               std::uint32_t bin_size, const std::string& assembly,
                               bool skip_all_vs_all_matrix, std::size_t threads,
                               std::size_t batch_size, std::uint32_t compression_lvl, bool force) {
  SPDLOG_INFO(FMT_STRING("begin loading pixels into a .hic file..."));

  if (force) {
    std::filesystem::remove(uri);
  }

  hic::internal::HiCFileWriter hf(uri, chromosomes, {bin_size}, assembly, threads, batch_size,
                                  tmp_dir, compression_lvl, skip_all_vs_all_matrix);

  std::vector<ThinPixel<float>> write_buffer(batch_size);
  return ingest_pixels(std::move(hf), pixel_queue, early_return, write_buffer);
}

inline Stats ingest_pairs_hic(PixelQueue<float>& pixel_queue, const std::atomic<bool>& early_return,
                              std::string_view uri, const std::filesystem::path& tmp_dir,
                              const Reference& chromosomes, std::uint32_t bin_size,
                              const std::string& assembly, bool skip_all_vs_all_matrix,
                              std::size_t threads, std::size_t batch_size,
                              std::uint32_t compression_lvl, bool force) {
  if (force) {
    std::filesystem::remove(uri);
  }

  hic::internal::HiCFileWriter hf(uri, chromosomes, {bin_size}, assembly, threads, batch_size,
                                  tmp_dir, compression_lvl, skip_all_vs_all_matrix);

  std::vector<ThinPixel<float>> buffer(batch_size);
  return ingest_pairs(std::move(hf), pixel_queue, early_return, buffer);
}

}  // namespace hictk::tools
