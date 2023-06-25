// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#include <algorithm>
#include <queue>
#include <string_view>
#include <vector>

#include "hictk/cooler.hpp"

namespace hictk::cooler::utils {

namespace internal {

[[nodiscard]] inline std::uint32_t get_bin_size_checked(const std::vector<File>& coolers) {
  assert(coolers.size() > 1);
  const auto& clr1 = coolers.front();

  for (std::size_t i = 1; i < coolers.size(); ++i) {
    const auto& clr2 = coolers[i];
    if (clr1.bin_size() != clr2.bin_size()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "cooler \"{}\" and \"{}\" have different resolutions ({}  and {} respectively)"),
          clr1.uri(), clr2.uri(), clr1.bin_size(), clr2.bin_size()));
    }
  }
  return clr1.bin_size();
}

[[nodiscard]] inline Reference get_chromosomes_checked(const std::vector<File>& coolers) {
  assert(coolers.size() > 1);
  const auto& clr1 = coolers.front();

  for (std::size_t i = 1; i < coolers.size(); ++i) {
    const auto& clr2 = coolers[i];
    if (clr1.chromosomes() != clr2.chromosomes()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("cooler \"{}\" and \"{}\" use different reference genomes"),
                      clr1.uri(), clr2.uri()));
    }
  }
  return clr1.chromosomes();
}

[[nodiscard]] inline bool merging_requires_float_pixels(const std::vector<File>& coolers) {
  for (const auto& clr : coolers) {
    if (clr.has_float_pixels()) {
      return true;
    }
  }
  return false;
}

template <typename N>
inline void merge(const std::vector<typename PixelSelector<N>::iterator>& heads,
                  const std::vector<typename PixelSelector<N>::iterator>& tails, File& dest,
                  std::size_t queue_capacity, bool quiet) {
  hictk::internal::PixelMerger merger{heads, tails};

  std::vector<Pixel<N>> buffer(queue_capacity);
  buffer.clear();

  std::size_t pixels_processed{};
  while (true) {
    auto pixel = merger.next();
    if (!pixel) {
      break;
    }

    buffer.emplace_back(std::move(pixel));
    if (buffer.size() == queue_capacity) {
      dest.append_pixels(buffer.begin(), buffer.end());
      pixels_processed += buffer.size();
      if (!quiet && pixels_processed % (std::max)(queue_capacity, std::size_t(10'000'000)) == 0) {
        spdlog::info(FMT_STRING("Procesed {}M pixels...\n"), pixels_processed / 10'000'000);
      }
    }
  }

  if (!buffer.empty()) {
    dest.append_pixels(buffer.begin(), buffer.end());
  }
}

template <typename N>
struct CoolerIteratorPairs {
  std::vector<typename PixelSelector<N>::iterator> heads{};
  std::vector<typename PixelSelector<N>::iterator> tails{};
};

template <typename N>
inline CoolerIteratorPairs<N> collect_iterators(const std::vector<File>& clrs) {
  if constexpr (std::is_floating_point_v<N>) {
    std::vector<PixelSelector<double>::iterator> heads{};
    std::vector<PixelSelector<double>::iterator> tails{};

    for (const auto& clr : clrs) {
      auto first = clr.begin<double>();
      auto last = clr.end<double>();
      if (first != last) {
        heads.emplace_back(std::move(first));
        tails.emplace_back(std::move(last));
      }
    }

    return {heads, tails};
  } else {
    std::vector<PixelSelector<std::int32_t>::iterator> heads{};
    std::vector<PixelSelector<std::int32_t>::iterator> tails{};

    for (const auto& clr : clrs) {
      auto first = clr.begin<std::int32_t>();
      auto last = clr.end<std::int32_t>();
      if (first != last) {
        heads.emplace_back(std::move(first));
        tails.emplace_back(std::move(last));
      }
    }
    return {heads, tails};
  }
}

}  // namespace internal

template <typename Str>
inline void merge(Str first_file, Str last_file, std::string_view dest_uri,
                  bool overwrite_if_exists, std::size_t chunk_size, bool quiet) {
  static_assert(std::is_constructible_v<std::string, decltype(*first_file)>);
  assert(chunk_size != 0);

  std::vector<File> clrs{};
  std::transform(first_file, last_file, std::back_inserter(clrs),
                 [&](const auto& uri) { return File::open_read_only_read_once(std::string{uri}); });

  if (clrs.size() < 2) {
    throw std::runtime_error("unable to merge less than 2 coolers");
  }

  const auto chroms = internal::get_chromosomes_checked(clrs);
  const auto bin_size = internal::get_bin_size_checked(clrs);
  const auto float_pixels = internal::merging_requires_float_pixels(clrs);

  auto dest =
      float_pixels
          ? File::create_new_cooler<double>(dest_uri, chroms, bin_size, overwrite_if_exists)
          : File::create_new_cooler<std::int32_t>(dest_uri, chroms, bin_size, overwrite_if_exists);
  try {
    if (float_pixels) {
      auto [heads, tails] = internal::collect_iterators<double>(clrs);
      internal::merge<double>(heads, tails, dest, chunk_size, quiet);
    } else {
      auto [heads, tails] = internal::collect_iterators<std::int32_t>(clrs);
      internal::merge<std::int32_t>(heads, tails, dest, chunk_size, quiet);
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to merge {} cooler files: {}"), clrs.size(), e.what()));
  }
}
}  // namespace hictk::cooler::utils
