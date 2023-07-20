// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/multires_cooler.hpp"
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

void zoomify_once(const cooler::File& clr1, cooler::RootGroup entrypoint2,
                  std::uint32_t resolution) {
  auto clr2 = cooler::File::create(std::move(entrypoint2), clr1.chromosomes(), resolution);

  std::vector<ThinPixel<std::int32_t>> buffer{500'000};
  cooler::MultiResFile::coarsen(clr1, clr2, buffer);
}

void zoomify_once(std::string_view uri1, std::string_view uri2, std::uint32_t resolution,
                  bool force) {
  auto clr1 = cooler::File::open(uri1);

  spdlog::info(FMT_STRING("coarsening cooler at {} once ({} -> {})"), clr1.uri(), clr1.bin_size(),
               resolution);

  auto mode = force ? HighFive::File::Overwrite : HighFive::File::Create;
  cooler::RootGroup entrypoint2{HighFive::File(std::string{uri2}, mode).getGroup("/")};

  return zoomify_once(clr1, std::move(entrypoint2), resolution);
}

void zoomify_many(std::string_view in_uri, std::string_view out_path,
                  const std::vector<std::uint32_t>& resolutions, bool copy_base_resolution,
                  bool force) {
  auto clr = cooler::File::open(in_uri);
  auto mclr =
      cooler::MultiResFile::create(out_path, cooler::File::open(in_uri).chromosomes(), force);

  spdlog::info(FMT_STRING("coarsening cooler at {} {} times ({} -> {})"), clr.uri(),
               resolutions.size(), clr.bin_size(), fmt::join(resolutions, " -> "));

  if (copy_base_resolution) {
    assert(resolutions.front() == clr.bin_size());
    mclr.copy_resolution(clr);
  } else {
    assert(resolutions.size() > 1);
    zoomify_once(cooler::File::open(in_uri), mclr.init_resolution(resolutions[1]), resolutions[1]);
  }

  for (std::size_t i = 1; i < resolutions.size(); ++i) {
    mclr.create_resolution(resolutions[i]);
  }
}

int zoomify_subcmd(const ZoomifyConfig& c) {
  const auto output_is_multires = c.copy_base_resolution || c.resolutions.size() > 2;

  if (output_is_multires) {
    zoomify_many(c.input_uri, c.output_path, c.resolutions, c.copy_base_resolution, c.force);
  } else {
    zoomify_once(c.input_uri, c.output_path, c.resolutions.back(), c.force);
  }

  return 0;
}
}  // namespace hictk::tools
