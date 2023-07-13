// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/transformers/coarsen.hpp"

namespace hictk::tools {

template <typename N = std::int32_t, typename PixelIt = cooler::PixelSelector<>::iterator<N>,
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

int zoomify_subcmd(const ZoomifyConfig& c) {
  cooler::init_mcool(c.output_path, c.resolutions.begin(), c.resolutions.end(), c.force);

  spdlog::info(FMT_STRING("Copying {} resolution from {}"), c.resolutions.front(), c.input_uri);
  cooler::utils::copy(c.input_uri, fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_path,
                                               c.resolutions.front()));

  std::vector<ThinPixel<std::int32_t>> buffer{500'000};

  const internal::TmpDir tmpdir{};
  for (std::size_t i = 1; i < c.resolutions.size(); ++i) {
    const auto res = c.resolutions[i];
    const auto base_resolution = [&]() {
      for (std::size_t j = i - 1; j != 0; --j) {
        if (res % c.resolutions[j] == 0) {
          return c.resolutions[j];
        }
      }
      return c.resolutions.front();
    }();

    const auto tmpcooler = tmpdir() / fmt::format(FMT_STRING("{}.cool"), res);
    const auto outcooler = fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_path, res);

    spdlog::info(FMT_STRING("Generating {} resolution from {} ({}x)"), res, base_resolution,
                 res / base_resolution);
    {
      auto clr1 = cooler::File::open_read_only_read_once(
          fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_path, base_resolution));
      auto clr2 = cooler::File::create_new_cooler(tmpcooler.string(), clr1.chromosomes(), res);

      auto merger = setup_pixel_merger(clr1, res / base_resolution);

      auto pixel = merger.next();
      buffer.clear();

      while (pixel) {
        buffer.emplace_back(std::move(pixel));
        if (buffer.size() == buffer.capacity()) {
          clr2.append_pixels(buffer.begin(), buffer.end());
          buffer.clear();
        }
        pixel = merger.next();
      }
      if (!buffer.empty()) {
        clr2.append_pixels(buffer.begin(), buffer.end());
      }
    }
    cooler::utils::copy(tmpcooler.string(), outcooler);
  }

  return 0;
}
}  // namespace hictk::tools
