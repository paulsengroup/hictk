// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/file_writer.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
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
  const auto path1 = (datadir / "4DNFIZ1ZVXC8.hic9").string();
  const auto path2 = (testdir() / "hic_block_partitioner.bin").string();
  const std::uint32_t resolution = 25'000;

  const hic::File f1(path1, resolution);
  const auto sel1 = f1.fetch("chr2L");
  const auto sel2 = f1.fetch("chr2L", "chr2R");

  const std::vector<ThinPixel<float>> pixels1(sel1.begin<float>(), sel1.end<float>());
  const std::vector<ThinPixel<float>> pixels2(sel2.begin<float>(), sel2.end<float>());

  HiCInteractionToBlockMapper partitioner(path2, f1.bins_ptr(), 50'000, 3);

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

TEST_CASE("HiC: HiCFileWriter", "[hic][v9][short]") {
  const auto path1 = (datadir / "4DNFIZ1ZVXC8.hic9").string();
  const auto path2 = (testdir() / "hic_writer_001.hic").string();
  const auto path3 = (testdir() / "hic_writer_002.hic").string();
  const std::vector<std::uint32_t> resolutions{500'000, 1'000'000, 2'500'000};

  SECTION("create file") {
    {
      const auto chromosomes = hic::File(path1, resolutions.front()).chromosomes();
      HiCFileWriter w(path2, chromosomes, resolutions, "dm6", 3);
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

  SECTION("add weights") {
    const std::uint32_t resolution = 500'000;
    const hic::File hf1(path1, resolution);

    {
      // init file
      HiCFileWriter w(path3, hf1.chromosomes(), {hf1.resolution()}, "dm6");
      const auto sel = hf1.fetch();
      w.add_pixels(resolution, sel.begin<float>(), sel.end<float>());
      w.serialize();
    }

    // add normalization weights
    {
      HiCFileWriter w(path3);
      for (const auto& chrom : w.chromosomes()) {
        if (chrom.is_all()) {
          continue;
        }
        w.add_norm_vector("SCALE", chrom, "BP", hf1.resolution(),
                          hf1.normalization("SCALE", chrom));
      }
      w.write_norm_vectors();
      CHECK_THROWS(w.add_norm_vector("VC", w.chromosomes().at("chr2L"), "BP", hf1.resolution(),
                                     std::vector<float>{1, 2, 3}));
    }

    // compare
    const hic::File hf2(path3, resolution);
    const auto pixels1 = hf1.fetch(balancing::Method::SCALE()).read_all<float>();
    const auto pixels2 = hf2.fetch(balancing::Method::SCALE()).read_all<float>();

    REQUIRE(pixels1.size() == pixels2.size());
    for (std::size_t i = 0; i < pixels1.size(); ++i) {
      CHECK(pixels1[i].coords == pixels2[i].coords);
      if (std::isnan(pixels1[i].count)) {
        CHECK(std::isnan(pixels2[i].count));
      } else {
        CHECK_THAT(pixels1[i].count, Catch::Matchers::WithinRel(pixels2[i].count));
      }
    }
  }
}

}  // namespace hictk::hic::test::file_writer
