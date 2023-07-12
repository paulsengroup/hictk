// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/tmpdir.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/transformers/coarsen.hpp"

namespace hictk::tools {

int zoomify_subcmd(const ZoomifyConfig& c) {
  cooler::init_mcool(c.output_path, c.resolutions.begin(), c.resolutions.end());

  cooler::utils::copy(c.input_uri, fmt::format(FMT_STRING("{}::/resolutions/{}"), c.output_path,
                                               c.resolutions.front()));

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

      auto sel1 = clr1.fetch();
      auto sel2 = transformers::CoarsenPixels(sel1.begin<std::int32_t>(), sel1.end<std::int32_t>(),
                                              res / base_resolution);
      clr2.append_pixels(sel2.begin(), sel2.end());
    }
    cooler::utils::copy(tmpcooler.string(), outcooler);
  }

  return 0;
}
}  // namespace hictk::tools
