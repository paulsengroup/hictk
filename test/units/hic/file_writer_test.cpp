// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/file_writer.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cstdint>
#include <filesystem>
#include <string>

#include "hictk/chromosome.hpp"
#include "hictk/hic.hpp"
#include "hictk/reference.hpp"
#include "hictk/transformers/join_genomic_coords.hpp"
#include "tmpdir.hpp"

using namespace hictk::hic;

namespace hictk::hic::test::file_writer {

using namespace hictk::hic::internal;

TEST_CASE("HiC: HiCInteractionToBlockMapper::BlockMapper", "[hic][v9][short]") {
  SECTION("intra") {
    // Test case based on blocks fetched by running an instrumented version of
    // hictk dump test/data/hic/4DNFIZ1ZVXC8.hic9 --resolution 10000 --range chr3R:0-50000
    {
      const HiCInteractionToBlockMapper::BlockMapperIntra mapper(803, 4);
      CHECK(mapper(0, 0) == 0);
      CHECK(mapper(0, 100) == 0);
      CHECK(mapper(802, 802) == 0);
      CHECK(mapper(803, 803) == 1);
      CHECK(mapper(1038, 2137) == 1);
      CHECK(mapper(235, 1376) == 5);
      CHECK(mapper(8, 3203) == 5);
    }
  }

  SECTION("inter") {
    // Test case based on blocks fetched by running an instrumented version of
    // hictk dump test/data/hic/4DNFIZ1ZVXC8.hic9 --resolution 10000 --range chr3L:0-50000 --range2
    // chr3R:0-10000000
    {
      const HiCInteractionToBlockMapper::BlockMapperInter mapper(803, 4);
      CHECK(mapper(0, 0) == 0);
      CHECK(mapper(0, 100) == 0);
      CHECK(mapper(802, 802) == 0);
      CHECK(mapper(7, 803) == 4);
      CHECK(mapper(795, 1605) == 4);
    }

    {
      const HiCInteractionToBlockMapper::BlockMapperInter mapper(101, 1);
      CHECK(mapper(0, 0) == 0);
      CHECK(mapper(0, 99) == 0);
      CHECK(mapper(99, 99) == 0);
    }
  }
}

TEST_CASE("HiC: HiCInteractionToBlockMapper", "[hic][v9][short]") {
  const std::uint32_t resolution = 25'000;
  const hic::File f1((datadir / "4DNFIZ1ZVXC8.hic9").string(), resolution);
  const auto sel1 = f1.fetch("chr2L");
  const auto sel2 = f1.fetch("chr2L", "chr2R");

  const std::vector<ThinPixel<float>> pixels1(sel1.begin<float>(), sel1.end<float>());
  const std::vector<ThinPixel<float>> pixels2(sel2.begin<float>(), sel2.end<float>());

  HiCInteractionToBlockMapper partitioner(testdir() / "hic_block_partitioner.bin", f1.bins_ptr(),
                                          50'000, 3);

  partitioner.append_pixels(pixels1.begin(), pixels1.end());
  partitioner.append_pixels(pixels2.begin(), pixels2.end());
  partitioner.finalize();

  std::size_t num_interactions = 0;
  for (const auto& [bid, _] : partitioner.block_index()) {
    const auto blk = partitioner.merge_blocks(bid);
    num_interactions += static_cast<std::size_t>(blk.nRecords);
  }

  CHECK(num_interactions == pixels1.size() + pixels2.size());
}

TEST_CASE("HiC: HiCFileWriter", "[hic][v9][long]") {
  const std::vector<std::uint32_t> resolutions{100'000, 200'000, 500'000, 1'000'000};

  const auto chromosomes =
      hic::File((datadir / "4DNFIZ1ZVXC8.hic9").string(), resolutions.front()).chromosomes();

  const auto path = testdir() / "hic_writer.hic";
  {
    HiCFileWriter w(path.string(), chromosomes, resolutions, "hg38", 16);
    for (std::size_t i = 0; i < resolutions.size(); ++i) {
      if (i % 2 == 0) {
        const auto resolution = resolutions[i];
        const hic::File f((datadir / "4DNFIZ1ZVXC8.hic9").string(), resolution);
        const auto sel = f.fetch();
        w.add_pixels(resolution, sel.begin<float>(), sel.end<float>());
      }
    }
    w.serialize();
  }

  for (const auto& resolution : resolutions) {
    fmt::print(FMT_STRING("Comparing {}...\n"), resolution);
    const hic::File f1((datadir / "4DNFIZ1ZVXC8.hic9").string(), resolution);
    const hic::File f2(path.string(), resolution);
    const auto expected_pixels = f1.fetch().read_all<float>();
    const auto pixels = f2.fetch().read_all<float>();

    REQUIRE(expected_pixels.size() == pixels.size());
    for (std::size_t i = 0; i < pixels.size(); ++i) {
      CHECK(expected_pixels[i] == pixels[i]);
    }
  }
}

}  // namespace hictk::hic::test::file_writer
