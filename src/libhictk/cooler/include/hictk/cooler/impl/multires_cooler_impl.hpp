// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <type_traits>
#include <variant>
#include <vector>

#include "hictk/cooler/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/pixel.hpp"
#include "hictk/transformers/coarsen.hpp"

namespace hictk::cooler {

template <typename ResolutionIt>
inline MultiResFile MultiResFile::create(const std::filesystem::path& path, const File& base,
                                         ResolutionIt first_res, ResolutionIt last_res,
                                         bool force_overwrite) {
  std::vector<std::uint32_t> resolutions_{first_res, last_res};
  std::sort(resolutions_.begin(), resolutions_.end());
  const auto base_res = base.resolution();
  const auto tgt_res = resolutions_.front();

  for (const auto& res : resolutions_) {
    if (base_res > tgt_res || res % base_res != 0) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("resolution {} is not a multiple of base resolution {}"), res, base_res));
    }
  }

  auto mclr = MultiResFile::create(path, base.chromosomes(), force_overwrite);
  SPDLOG_INFO(FMT_STRING("Copying {} resolution from \"{}\""), base_res, base.path());
  mclr.copy_resolution(base);

  std::visit(
      [&]([[maybe_unused]] auto count_type) {
        using PixelT = decltype(count_type);

        std::vector<ThinPixel<PixelT>> buffer{500'000};  // NOLINT(*-avoid-magic-numbers)
        for (std::size_t i = 1; i < resolutions_.size(); ++i) {
          const auto tgt_resolution = resolutions_[i];
          const auto base_resolution = compute_base_resolution(mclr._resolutions, tgt_resolution);

          auto attributes = Attributes::init<PixelT>(tgt_resolution);
          attributes.assembly = base.attributes().assembly;

          auto clr = File::create<PixelT>(mclr.init_resolution(tgt_resolution), base.chromosomes(),
                                          tgt_resolution, attributes, DEFAULT_HDF5_CACHE_SIZE);

          SPDLOG_INFO(FMT_STRING("Generating {} resolution from {} ({}x)"), tgt_resolution,
                      base_resolution, tgt_resolution / base_resolution);

          MultiResFile::coarsen(mclr.open(base_res), clr, buffer);
          mclr._resolutions.push_back(tgt_resolution);
        }
      },
      base.pixel_variant());

  return mclr;
}

constexpr const std::vector<std::uint32_t>& MultiResFile::resolutions() const noexcept {
  return _resolutions;
}

constexpr const MultiResAttributes& MultiResFile::attributes() const noexcept { return _attrs; }

template <typename N>
inline File MultiResFile::create_resolution(std::uint32_t resolution, Attributes attributes) {
  const auto base_resolution = compute_base_resolution(resolutions(), resolution);

  std::vector<ThinPixel<N>> buffer{500'000};  // NOLINT(*-avoid-magic-numbers)
  auto base_clr = open(base_resolution);
  attributes.assembly = base_clr.attributes().assembly;
  attributes.bin_size = resolution;
  {
    auto clr = File::create<N>(init_resolution(resolution), base_clr.chromosomes(), resolution,
                               attributes, DEFAULT_HDF5_CACHE_SIZE);
    MultiResFile::coarsen(base_clr, clr, buffer);
  }

  _resolutions.push_back(resolution);
  std::sort(_resolutions.begin(), _resolutions.end());

  return open(resolution);
}

template <typename N>
inline void MultiResFile::coarsen(const File& clr1, File& clr2, std::vector<ThinPixel<N>>& buffer) {
  static_assert(std::is_arithmetic_v<N>);
  SPDLOG_INFO(FMT_STRING("generating {} resolution from {} ({}x)"), clr2.resolution(),
              clr1.resolution(), clr2.resolution() / clr1.resolution());
  auto sel1 = clr1.fetch();
  auto sel2 = transformers::CoarsenPixels(sel1.begin<N>(), sel1.end<N>(), clr1.bins_ptr(),
                                          clr2.resolution() / clr1.resolution());

  const auto update_frequency =
      std::max(std::size_t{1'000'000}, (clr1.dataset("pixels/bin1_id").size() / 100));

  auto first = sel2.begin();
  auto last = sel2.end();
  buffer.clear();

  auto t0 = std::chrono::steady_clock::now();
  for (std::size_t j = 0; first != last; ++j) {
    buffer.emplace_back(*first);
    if (buffer.size() == buffer.capacity()) {
      clr2.append_pixels(buffer.begin(), buffer.end());
      buffer.clear();
    }
    if (j == update_frequency) {
      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      const auto bin1 = clr2.bins().at(first->bin1_id);
      SPDLOG_INFO(FMT_STRING("[{} -> {}] processing {:ucsc} at {:.0f} pixels/s..."),
                  clr1.resolution(), clr2.resolution(), bin1,
                  static_cast<double>(update_frequency) / delta);
      t0 = t1;
      j = 0;
    }
    ++first;
  }
  if (!buffer.empty()) {
    clr2.append_pixels(buffer.begin(), buffer.end());
  }
}

}  // namespace hictk::cooler
