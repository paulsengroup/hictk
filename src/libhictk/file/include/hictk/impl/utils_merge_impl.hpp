// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/transformers/pixel_merger.hpp"

namespace hictk::utils {

namespace merge::internal {

inline void validate_bin_size(const std::vector<File>& files, bool variable_bin_sizes_ok = true) {
  assert(files.size() > 1);
  const auto& f1 = files.front();
  const auto bin_table1 = f1.bins();
  if (!variable_bin_sizes_ok && f1.bins().type() == BinTable::Type::variable) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("file \"{}\" has a bin table with variable bin size"), f1.uri()));
  }

  for (std::size_t i = 1; i < files.size(); ++i) {
    const auto& f2 = files[i];

    if (f1.bins().type() == BinTable::Type::variable ||
        f2.bins().type() == BinTable::Type::variable) {
      if (!variable_bin_sizes_ok && f2.bins().type() == BinTable::Type::variable) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("file \"{}\" has a bin table with variable bin size"), f2.uri()));
      }
      const auto& bin_table2 = f2.bins();
      if (bin_table1 != bin_table2) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("files \"{}\" and \"{}\" have different bin tables"), f1.uri(), f2.uri()));
      }
    } else if (f1.resolution() != f2.resolution()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("files \"{}\" and \"{}\" have different resolutions ({} and {} respectively)"),
          f1.uri(), f2.uri(), f1.resolution(), f2.resolution()));
    }
  }
}

inline void validate_chromosomes(const std::vector<File>& files) {
  assert(files.size() > 1);
  const auto& f1 = files.front();

  for (std::size_t i = 1; i < files.size(); ++i) {
    const auto& f2 = files[i];
    if (f1.chromosomes() != f2.chromosomes()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("files \"{}\" and \"{}\" use different reference genomes"),
                      f1.uri(), f2.uri()));
    }
  }
}

template <typename N>
[[nodiscard]] inline std::tuple<std::vector<PixelSelector>, std::vector<PixelSelector::iterator<N>>,
                                std::vector<PixelSelector::iterator<N>>>
init_iterators(const std::vector<File>& files) {
  std::vector<PixelSelector> selectors{};
  std::vector<PixelSelector::iterator<N>> heads{};
  std::vector<PixelSelector::iterator<N>> tails{};
  for (const auto& f : files) {
    const auto& sel = selectors.emplace_back(f.fetch());
    auto first = sel.begin<N>();
    auto last = sel.end<N>();
    if (first == last) {
      selectors.pop_back();
    } else {
      heads.emplace_back(std::move(first));
      tails.emplace_back(std::move(last));
    }
  }

  return std::make_tuple(std::move(selectors), std::move(heads), std::move(tails));
}

}  // namespace merge::internal

/// Iterable of strings
template <typename N, typename Str>
inline void merge_to_cool(Str first_uri, Str last_uri, std::string_view dest_uri,
                          std::uint32_t resolution, bool overwrite_if_exists,
                          std::size_t chunk_size, std::size_t update_frequency,
                          std::uint32_t compression_lvl) {
  if (std::distance(first_uri, last_uri) < 2) {
    throw std::runtime_error("cannot merge less than 2 files");
  }
  if (std::all_of(first_uri, last_uri,
                  [&](const auto& uri) { return cooler::utils::is_cooler(uri); })) {
    return cooler::utils::merge<N>(first_uri, last_uri, dest_uri, overwrite_if_exists, chunk_size,
                                   update_frequency, compression_lvl);
  }

  std::vector<File> files{};
  std::transform(first_uri, last_uri, std::back_inserter(files),
                 [&](const auto& uri) { return File(uri, resolution); });
  try {
    merge::internal::validate_chromosomes(files);
    merge::internal::validate_bin_size(files);
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("cannot merge files {}: {}"),
                                         fmt::join(first_uri, last_uri, ", "), e.what()));
  }

  const auto& [_, heads, tails] = merge::internal::init_iterators<N>(files);

  transformers::PixelMerger merger{heads, tails};
  std::vector<ThinPixel<N>> buffer(chunk_size);
  buffer.clear();

  const auto& f0 = files.front();
  auto attrs = cooler::Attributes::init(f0.resolution());
  std::visit(
      [&](const auto& f) {
        using FileT = remove_cvref_t<decltype(f)>;
        if constexpr (std::is_same_v<FileT, cooler::File>) {
          attrs.assembly = f.attributes().assembly;
        } else {
          attrs.assembly = f.assembly();
        }
      },
      f0.get());

  auto dest = cooler::File::create<N>(dest_uri, f0.bins(), overwrite_if_exists, attrs,
                                      cooler::DEFAULT_HDF5_CACHE_SIZE * 4, compression_lvl);

  dest.append_pixels(merger.begin(), merger.end());
}

/// Iterable of strings
template <typename Str>
inline void merge_to_hic(Str first_file, Str last_file, std::string_view dest_file,
                         std::uint32_t resolution, const std::filesystem::path& tmp_dir,
                         bool overwrite_if_exists, std::size_t chunk_size, std::size_t n_threads,
                         std::uint32_t compression_lvl, bool skip_all_vs_all) {
  if (std::distance(first_file, last_file) < 2) {
    throw std::runtime_error("cannot merge less than 2 files");
  }
  if (std::all_of(first_file, last_file,
                  [&](const auto& uri) { return hic::utils::is_hic_file(uri); })) {
    return hic::utils::merge(first_file, last_file, dest_file, resolution, tmp_dir,
                             overwrite_if_exists, chunk_size, n_threads, compression_lvl,
                             skip_all_vs_all);
  }

  std::vector<File> files{};
  std::transform(first_file, last_file, std::back_inserter(files),
                 [&](const auto& uri) { return File(uri, resolution); });
  try {
    merge::internal::validate_chromosomes(files);
    merge::internal::validate_bin_size(files, false);
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("cannot merge files {}: {}"),
                                         fmt::join(first_file, last_file, ", "), e.what()));
  }

  if (overwrite_if_exists) {
    std::filesystem::remove(dest_file);
  }

  const auto& [_, heads, tails] = merge::internal::init_iterators<float>(files);

  transformers::PixelMerger merger{heads, tails};
  std::vector<ThinPixel<float>> buffer(chunk_size);
  buffer.clear();

  const auto& f0 = files.front();
  const auto assembly = std::visit(
      [&](const auto& f) -> std::string {
        using FileT = remove_cvref_t<decltype(f)>;
        if constexpr (std::is_same_v<FileT, cooler::File>) {
          return std::string{f.attributes().assembly.has_value() ? *f.attributes().assembly
                                                                 : "unknown"};
        } else {
          return std::string{f.assembly()};
        }
      },
      f0.get());

  hic::internal::HiCFileWriter w(dest_file, f0.chromosomes(), {f0.resolution()}, assembly,
                                 n_threads, chunk_size, tmp_dir, compression_lvl, skip_all_vs_all);

  w.add_pixels(f0.resolution(), merger.begin(), merger.end());
  w.serialize();
}

}  // namespace hictk::utils
