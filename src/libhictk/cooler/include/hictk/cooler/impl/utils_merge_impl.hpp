// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#include <algorithm>
#include <queue>
#include <string_view>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/transformers/pixel_merger.hpp"

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
  auto clr = File::open_read_once(std::string{uri});
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

}  // namespace internal

template <typename N, typename Str>
inline void merge(Str first_uri, Str last_uri, std::string_view dest_uri, bool overwrite_if_exists,
                  std::size_t chunk_size, std::size_t update_frequency) {
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

    std::vector<PixelSelector::iterator<N>> heads;
    std::vector<PixelSelector::iterator<N>> tails;

    for (auto& clr : clrs) {
      if (clr.first_pixel != clr.last_pixel) {
        heads.emplace_back(std::move(clr.first_pixel));
        tails.emplace_back(std::move(clr.last_pixel));
      }
    }

    merge(heads, tails, cooler::File(clrs.front().uri).bins(), dest_uri, overwrite_if_exists,
          chunk_size, update_frequency);
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to merge {} cooler files: {}"),
                                         std::distance(first_uri, last_uri), e.what()));
  }
}

template <typename PixelIt>
inline void merge(const std::vector<PixelIt>& heads, const std::vector<PixelIt>& tails,
                  const BinTable& bins, std::string_view dest_uri, bool overwrite_if_exists,
                  std::size_t chunk_size, std::size_t update_frequency) {
  using N = remove_cvref_t<decltype(heads.front()->count)>;

  hictk::transformers::PixelMerger merger{heads, tails};
  std::vector<ThinPixel<N>> buffer(chunk_size);
  buffer.clear();

  auto dest = File::create<N>(dest_uri, bins, overwrite_if_exists);

  std::size_t pixels_processed{};
  auto t0 = std::chrono::steady_clock::now();
  auto first = merger.begin();
  auto last = merger.end();
  for (std::size_t i = 0; first != last; ++i) {
    const auto pixel = *first;
    std::ignore = ++first;

    if (i == update_frequency) {
      const auto bin1 = dest.bins().at(pixel.bin1_id);
      const auto bin2 = dest.bins().at(pixel.bin2_id);

      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      SPDLOG_INFO(FMT_STRING("processing {:ucsc} {:ucsc} at {:.0f} pixels/s..."), bin1, bin2,
                  double(update_frequency) / delta);
      t0 = t1;
      i = 0;
    }

    buffer.emplace_back(std::move(pixel));
    if (buffer.size() == chunk_size) {
      dest.append_pixels(buffer.begin(), buffer.end());
      pixels_processed += buffer.size();
      buffer.clear();
    }
  }

  if (!buffer.empty()) {
    dest.append_pixels(buffer.begin(), buffer.end());
  }
}

}  // namespace hictk::cooler::utils
