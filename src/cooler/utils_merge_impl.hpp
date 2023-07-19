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

template <typename N>
struct LightCooler {
  std::string uri{};
  Reference chroms{};
  std::uint32_t bin_size{};

  PixelSelector::iterator<N> first_pixel{};
  PixelSelector::iterator<N> last_pixel{};
};

template <typename N>
[[nodiscard]] inline LightCooler<N> preprocess_cooler(const std::string& uri) {
  auto clr = File::open_read_only_read_once(std::string{uri});
  auto sel = clr.fetch();
  return {uri, clr.chromosomes(), clr.bin_size(), sel.begin<N>(), sel.end<N>()};
}

template <typename N>
inline void validate_bin_size(const std::vector<LightCooler<N>>& coolers) {
  assert(coolers.size() > 1);
  const auto& clr1 = coolers.front();

  for (std::size_t i = 1; i < coolers.size(); ++i) {
    const auto& clr2 = coolers[i];
    if (clr1.bin_size != clr2.bin_size) {
      throw std::runtime_error(fmt::format(
          FMT_STRING(
              "cooler \"{}\" and \"{}\" have different resolutions ({}  and {} respectively)"),
          clr1.uri, clr2.uri, clr1.bin_size, clr2.bin_size));
    }
  }
}

template <typename N>
inline void validate_chromosomes(const std::vector<LightCooler<N>>& coolers) {
  assert(coolers.size() > 1);
  const auto& clr1 = coolers.front();

  for (std::size_t i = 1; i < coolers.size(); ++i) {
    const auto& clr2 = coolers[i];
    if (clr1.chroms != clr2.chroms) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("cooler \"{}\" and \"{}\" use different reference genomes"),
                      clr1.uri, clr2.uri));
    }
  }
}

template <typename N>
inline void merge(const std::vector<PixelSelector::iterator<N>>& heads,
                  const std::vector<PixelSelector::iterator<N>>& tails, File& dest,
                  std::size_t queue_capacity, bool quiet) {
  hictk::internal::PixelMerger merger{heads, tails};

  std::vector<ThinPixel<N>> buffer(queue_capacity);
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
      buffer.clear();
    }
  }

  if (!buffer.empty()) {
    dest.append_pixels(buffer.begin(), buffer.end());
  }
}

}  // namespace internal

template <typename N, typename Str>
inline void merge(Str first_uri, Str last_uri, std::string_view dest_uri, bool overwrite_if_exists,
                  std::size_t chunk_size, bool quiet) {
  static_assert(std::is_constructible_v<std::string, decltype(*first_uri)>);
  assert(chunk_size != 0);
  try {
    std::vector<internal::LightCooler<N>> clrs{};
    std::transform(first_uri, last_uri, std::back_inserter(clrs), [&](const auto& uri) {
      return internal::preprocess_cooler<N>(std::string{uri});
    });

    if (clrs.size() < 2) {
      throw std::runtime_error("cannot merge less than 2 coolers");
    }

    internal::validate_chromosomes(clrs);
    internal::validate_bin_size(clrs);
    auto dest = File::create_new_cooler<N>(dest_uri, clrs.front().chroms, clrs.front().bin_size,
                                           overwrite_if_exists);

    std::vector<PixelSelector::iterator<N>> heads;
    std::vector<PixelSelector::iterator<N>> tails;

    for (auto& clr : clrs) {
      if (clr.first_pixel != clr.last_pixel) {
        heads.emplace_back(std::move(clr.first_pixel));
        tails.emplace_back(std::move(clr.last_pixel));
      }
    }

    internal::merge(heads, tails, dest, chunk_size, quiet);
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to merge {} cooler files: {}"),
                                         std::distance(first_uri, last_uri), e.what()));
  }
}
}  // namespace hictk::cooler::utils
