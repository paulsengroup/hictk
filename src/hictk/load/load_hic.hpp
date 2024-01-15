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
                               const std::string& assembly, std::int64_t offset, Format format,
                               std::size_t threads, std::size_t batch_size,
                               std::int32_t compression_lvl, bool force) {
  SPDLOG_INFO(FMT_STRING("begin loading pixels into a .hic file..."));

  if (force) {
    std::filesystem::remove(uri);
  }

  hic::internal::HiCHeader header{
      std::string{uri},  // url
      9,                 // version
      -1,                // masterIndexOffset
      assembly,          // genomeID
      -1,                // nviPosition
      -1,                // nviLength
      chromosomes,
      {bin_size},              // resolutions
      {{"software", "hictk"}}  // attributes
  };

  hic::internal::HiCFileWriter hf(std::move(header), threads, batch_size, tmp_dir, compression_lvl);

  std::vector<ThinPixel<float>> write_buffer(batch_size);
  return ingest_pixels(std::move(hf), write_buffer, format, offset);
}

inline Stats ingest_pairs_hic(std::string_view uri, const std::filesystem::path& tmp_dir,
                              const Reference& chromosomes, std::uint32_t bin_size,
                              const std::string& assembly, std::int64_t offset, Format format,
                              std::size_t threads, std::size_t batch_size,
                              std::int32_t compression_lvl, bool force) {
  if (force) {
    std::filesystem::remove(uri);
  }

  hic::internal::HiCHeader header{
      std::string{uri},  // url
      9,                 // version
      -1,                // masterIndexOffset
      assembly,          // genomeID
      -1,                // nviPosition
      -1,                // nviLength
      chromosomes,
      {bin_size},              // resolutions
      {{"software", "hictk"}}  // attributes
  };

  hic::internal::HiCFileWriter hf(std::move(header), threads, batch_size, tmp_dir, compression_lvl);

  std::vector<ThinPixel<float>> buffer(batch_size);
  return ingest_pairs(std::move(hf), buffer, format, offset);
}

}  // namespace hictk::tools
