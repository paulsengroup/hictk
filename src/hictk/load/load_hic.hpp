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

  const auto resolution = hf.resolutions().front();
  assert(buffer.capacity() != 0);
  buffer.reserve(buffer.capacity());

  for (std::size_t i = 0; true; ++i) {
    SPDLOG_INFO(FMT_STRING("preprocessing chunk #{}..."), i + 1);
    PairsAggregator<float>{hf.bins(resolution), format, offset}.read_next_chunk(buffer);

    if (buffer.empty()) {
      assert(std::cin.eof());
      break;
    }

    hf.add_pixels(resolution, buffer.begin(), buffer.end());
    SPDLOG_INFO(FMT_STRING("done preprocessing chunk #{}"), i + 1);
  }

  hf.serialize();
  const auto stats = hf.stats(resolution);
  return {stats.sum, stats.nnz};
}

}  // namespace hictk::tools
