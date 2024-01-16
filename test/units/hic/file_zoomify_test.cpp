// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/file_zoomify.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "tmpdir.hpp"

using namespace hictk::hic;

namespace hictk::hic::test::file_writer {

using namespace hictk::hic::internal;

TEST_CASE("HiC: HiCFileZoomify") {
  const std::vector<std::uint32_t> resolutions{100'000, 400'000, 1'000'000};

  const auto path1 = (datadir / "4DNFIZ1ZVXC8.hic9").string();
  const auto path2 = (testdir() / "hic_file_zoomify.hic").string();

  const auto avail_resolutions = hic::utils::list_resolutions(path1);
  REQUIRE(std::find(avail_resolutions.begin(), avail_resolutions.end(), 400'000) ==
          avail_resolutions.end());

  {
    hic::internal::HiCFileZoomify hzmf(path1, path2, resolutions);
    hzmf.zoomify();
  }
  for (const auto& resolution : {100'000U, 1'000'000U}) {
    const hic::File f1(path1, resolution);
    const hic::File f2(path2, resolution);

    const auto expected_pixels = f1.fetch().read_all<float>();
    const auto pixels = f2.fetch().read_all<float>();

    REQUIRE(expected_pixels.size() == pixels.size());
    for (std::size_t i = 0; i < pixels.size(); ++i) {
      CHECK(expected_pixels[i] == pixels[i]);
    }
  }
}

}  // namespace hictk::hic::test::file_writer
