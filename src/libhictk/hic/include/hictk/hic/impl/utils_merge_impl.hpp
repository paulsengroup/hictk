// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/file_writer.hpp"

namespace hictk::hic::utils {

namespace internal {

inline void validate_chromosomes(const std::vector<hic::File>& files) {
  assert(files.size() > 1);
  const auto chromosomes = files.front().chromosomes().remove_ALL();

  for (std::size_t i = 1; i < files.size(); ++i) {
    if (chromosomes != files[i].chromosomes().remove_ALL()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("files \"{}\" and \"{}\" use different reference genomes"),
                      files.front().path(), files[i].path()));
    }
  }
}
}  // namespace internal

template <typename Str>
inline void merge(Str first_file, Str last_file, std::string_view dest_file,
                  std::uint32_t resolution, const std::filesystem::path& tmp_dir,
                  bool overwrite_if_exists, std::size_t chunk_size, std::size_t n_threads,
                  std::uint32_t compression_lvl, bool skip_all_vs_all) {
  static_assert(std::is_constructible_v<std::string, decltype(*first_file)>);
  assert(chunk_size != 0);
  try {
    std::vector<hic::File> files{};
    std::transform(first_file, last_file, std::back_inserter(files),
                   [&](const std::string& path) { return hic::File(path, resolution); });
    if (files.size() < 2) {
      throw std::runtime_error("cannot merge less than 2 .hic files");
    }

    internal::validate_chromosomes(files);

    std::vector<PixelSelectorAll> selectors;
    std::vector<PixelSelectorAll::iterator<float>> heads;
    std::vector<PixelSelectorAll::iterator<float>> tails;

    for (auto& f : files) {
      auto sel = f.fetch();
      auto first = sel.begin<float>();
      auto last = sel.end<float>();
      if (first != last) {
        selectors.emplace_back(std::move(sel));
        heads.emplace_back(std::move(first));
        tails.emplace_back(std::move(last));
      }
    }

    merge(heads, tails, files.front().bins(), dest_file, files.front().assembly(), tmp_dir,
          overwrite_if_exists, chunk_size, n_threads, compression_lvl, skip_all_vs_all);
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to merge {} .hic files: {}"),
                                         std::distance(first_file, last_file), e.what()));
  }
}

template <typename PixelIt>
inline void merge(const std::vector<PixelIt>& heads, const std::vector<PixelIt>& tails,
                  const BinTable& bins, std::string_view dest_file, std::string_view assembly,
                  const std::filesystem::path& tmp_dir, bool overwrite_if_exists,
                  std::size_t chunk_size, std::size_t n_threads, std::uint32_t compression_lvl,
                  bool skip_all_vs_all) {
  using N = remove_cvref_t<decltype(heads.front()->count)>;

  const transformers::PixelMerger merger{heads, tails};
  std::vector<ThinPixel<N>> buffer(chunk_size);
  buffer.clear();

  if (overwrite_if_exists) {
    std::filesystem::remove(dest_file);
  }

  hic::internal::HiCFileWriter w(dest_file, bins.chromosomes(), {bins.resolution()}, assembly,
                                 n_threads, chunk_size, tmp_dir, compression_lvl, skip_all_vs_all);

  w.add_pixels(bins.resolution(), merger.begin(), merger.end());
  w.serialize();
}

}  // namespace hictk::hic::utils
