// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/transformers/coarsen.hpp"

namespace hictk::tools {

template <typename N = std::int32_t, typename PixelIt = cooler::PixelSelector::iterator<N>,
          typename CoarsenIt = typename transformers::CoarsenPixels<PixelIt>::iterator>
internal::PixelMerger<CoarsenIt> setup_pixel_merger(const cooler::File& clr, std::size_t factor) {
  const auto& chroms = clr.chromosomes();
  std::vector<CoarsenIt> heads{};
  std::vector<CoarsenIt> tails{};

  for (std::uint32_t chrom1_id = 0; chrom1_id < chroms.size(); ++chrom1_id) {
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chroms.size(); ++chrom2_id) {
      auto sel = clr.fetch(chroms.at(chrom1_id).name(), chroms.at(chrom2_id).name());
      auto sel1 = transformers::CoarsenPixels(sel.begin<N>(), sel.end<N>(), clr.bins_ptr(), factor);
      auto first = sel1.begin();
      auto last = sel1.end();
      if (first != last) {
        heads.emplace_back(std::move(first));
        tails.emplace_back(std::move(last));
      }
    }
  }
  return {heads, tails};
}

[[nodiscard]] static std::uint32_t compute_base_resolution(
    const std::vector<std::uint32_t>& resolutions, std::uint32_t target_res) {
  assert(!resolutions.empty());
  assert(target_res < resolutions.front());
  assert(target_res % resolutions.front() == 0);

  return *std::find_if(resolutions.rbegin(), resolutions.rend(),
                       [&](const auto res) { return res < target_res && res % target_res == 0; });
}

static void coarsen_cooler(std::string input_uri, std::string output_uri, std::uint32_t resolution,
                           std::vector<ThinPixel<std::int32_t>>& buffer) {
  auto clr1 = cooler::File::open_read_once(input_uri);
  auto clr2 = cooler::File::create(output_uri, clr1.chromosomes(), resolution);

  auto sel1 = clr1.fetch();
  auto sel2 = transformers::CoarsenPixels(sel1.begin<std::int32_t>(), sel1.end<std::int32_t>(),
                                          clr1.bins_ptr(), clr2.bin_size() / clr1.bin_size());

  const auto update_frequency =
      std::max(std::size_t(1'000'000), (clr1.dataset("pixels/bin1_id").size() / 100));

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
      spdlog::info(FMT_STRING("[{} -> {}] processing {:ucsc} at {:.0f} pixels/s..."),
                   clr1.bin_size(), resolution, bin1, double(update_frequency) / delta);
      t0 = t1;
      j = 0;
    }
    ++first;
  }
  if (!buffer.empty()) {
    clr2.append_pixels(buffer.begin(), buffer.end());
  }
}

int zoomify_subcmd(const ZoomifyConfig& c) {
  const auto output_is_multires = c.copy_base_resolution || c.resolutions.size() > 2;

  if (output_is_multires) {
    auto resolutions = c.resolutions;
    if (!c.copy_base_resolution) {
      resolutions.erase(resolutions.begin());
    }
    cooler::init_mcool(c.output_path, resolutions.begin(), resolutions.end(), c.force);
  }

  if (c.copy_base_resolution) {
    spdlog::info(FMT_STRING("Copying {} resolution from {}"), c.resolutions.front(), c.input_uri);
    cooler::utils::copy(c.input_uri, fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_path,
                                                 c.resolutions.front()));
  }

  const auto input_resolution = c.resolutions.front();

  std::vector<ThinPixel<std::int32_t>> buffer{500'000};

  const internal::TmpDir tmpdir{};  // TODO update
  for (std::size_t i = 1; i < c.resolutions.size(); ++i) {
    const auto tgt_resolution = c.resolutions[i];
    const auto base_resolution = compute_base_resolution(c.resolutions, tgt_resolution);

    const auto tmpcooler = tmpdir() / fmt::format(FMT_STRING("{}.cool"), tgt_resolution);
    const auto outcooler =
        fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_path, tgt_resolution);
    const auto incooler =
        input_resolution == base_resolution
            ? c.input_uri
            : fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_path, base_resolution);

    spdlog::info(FMT_STRING("Generating {} resolution from {} ({}x)"), tgt_resolution,
                 base_resolution, tgt_resolution / base_resolution);

    coarsen_cooler(incooler, tmpcooler.string(), tgt_resolution, buffer);
    cooler::utils::copy(tmpcooler.string(), outcooler);
    std::filesystem::remove(tmpcooler);
  }

  return 0;
}
}  // namespace hictk::tools
