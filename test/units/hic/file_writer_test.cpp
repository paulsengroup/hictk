// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef SPDLOG_ACTIVE_LEVEL
#undef SPDLOG_ACTIVE_LEVEL
#endif
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

#include "hictk/hic/file_writer.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <cstdint>
#include <filesystem>
#include <string>

#include "hictk/chromosome.hpp"
#include "hictk/fmt.hpp"  // TODO remove
#include "hictk/hic.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/reference.hpp"
#include "tmpdir.hpp"

using namespace hictk::hic;

namespace hictk::hic::test::file_writer {

using namespace hictk::hic::internal;

/*
TEST_CASE("HiC: write header", "[hic][v9][short]") {
  const std::int64_t master_index_offset = 115;

  // clang-format off
  const HiCHeader header1{
      testdir() / "hic.header",  // url
      9,                         // version
      master_index_offset,
      "hg38",                    // genomeID
      -1,                        // nviPosition
      -1,                        // nviLength
      Reference{Chromosome{0, "chr1", 248956422},
                Chromosome{1, "chr2", 242193529},
                Chromosome{2, "chr3", 198295559}},
      {100U, 1000U, 10'000},     // resolutions
      {{"software", "hictk"}}    // attributes
  };
  // clang-format on

  {
    HiCFileWriter w(header1);

    w.write_header();
    w.write_footer_offset(master_index_offset);

    std::ofstream ofs(header1.url, std::ios::binary | std::ios::app);
    ofs.write("0", 1);  // Add padding
  }

  const auto header2 = HiCFileReader{header1.url}.header();
  CHECK(header1 == header2);
}
*/

/*
TEST_CASE("HiC: BlockMapperIntra", "[hic][v9][short]") {
  SECTION("intra") {
    // Test case based on blocks fetched by running an instrumented version of
    // hictk dump test/data/hic/4DNFIZ1ZVXC8.hic9 --resolution 10000 --range chr3R:0-50000
    {
      const ChromChromHiCFileWriter::BlockMapperIntra mapper(803, 4);
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
      const ChromChromHiCFileWriter::BlockMapperInter mapper(803, 4);
      CHECK(mapper(0, 0) == 0);
      CHECK(mapper(0, 100) == 0);
      CHECK(mapper(802, 802) == 0);
      CHECK(mapper(7, 803) == 4);
      CHECK(mapper(795, 1605) == 4);
    }

    {
      const ChromChromHiCFileWriter::BlockMapperInter mapper(101, 1);
      CHECK(mapper(0, 0) == 0);
      CHECK(mapper(0, 99) == 0);
      CHECK(mapper(99, 99) == 0);
    }
  }
}
*/

/*
TEST_CASE("devel") {
  {
    // clang-format off
  const HiCHeader header{
          "/tmp/test.hic",           // url
          9,                         // version
          -1,                        // masterIndexOffset
          "hg38",                    // genomeID
          -1,                        // nviPosition
          -1,                        // nviLength
          Reference{Chromosome{0, "ALL", 123},
                    Chromosome{1, "chr1", 248956422},
                    Chromosome{2, "chr2", 242193529},
                    Chromosome{3, "chr3", 198295559}},
          {1000U},                   // resolutions
          {{"software", "hictk"}}    // attributes
  };
    // clang-format on

    HiCFileWriter w(header);
    w.write_header();

    std::vector<ThinPixel<float>> pixels;
    for (std::size_t i = 0; i < 100; ++i) {
      for (std::size_t j = i; j < 100; ++j) {
        pixels.push_back({i, j, 10});
      }
    }
    std::vector<HiCFooter> footers{};
    std::vector<std::int64_t> master_index_offsets{};
    std::vector<std::int32_t> matrix_metadata_offsets{};

    for (const auto& chrom : header.chromosomes) {
      const HiCFooter footer{
          {},
          HiCFooterMetadata{"/tmp/test.hic", MatrixType::observed, balancing::Method::NONE(),
                            MatrixUnit::BP, 1000, chrom, chrom, 0},
          {},
          {},
          {}

      };

      if (chrom == "chr1") {
        w.write_interaction_block(0, w.chromosomes().at("chr1"), w.chromosomes().at("chr1"), 1000,
                                  pixels, 0, 0);
      }

      const auto [master_index_offset, matrix_metadata_size] = w.write_body_metadata(chrom, chrom);
      master_index_offsets.push_back(master_index_offset);
      matrix_metadata_offsets.push_back(static_cast<std::int32_t>(matrix_metadata_size));
      footers.push_back(footer);
    }
    const auto footer_offset =
        w.write_footer(footers, master_index_offsets, matrix_metadata_offsets);
    w.write_master_index_offset(footer_offset);
  }

  hic::File f("/tmp/test.hic", 1000);
  auto sel = f.fetch("chr1");

  std::for_each(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), [](const auto& p) {
    fmt::print(FMT_STRING("{}\t{}\t{}\n"), p.bin1_id, p.bin2_id, p.count);
  });
}

*/

TEST_CASE("HiC: HiCInteractionToBlockMapper") {
  spdlog::default_logger()->set_level(spdlog::level::debug);
  const std::uint32_t resolution = 25'000;
  const hic::File f1((datadir / "4DNFIZ1ZVXC8.hic9").string(), resolution);
  const auto sel1 = f1.fetch("chr2L");
  const auto sel2 = f1.fetch("chr2L", "chr2R");

  const std::vector<ThinPixel<float>> pixels1(sel1.begin<float>(), sel1.end<float>());
  const std::vector<ThinPixel<float>> pixels2(sel2.begin<float>(), sel2.end<float>());

  HiCInteractionToBlockMapper partitioner(testdir() / "hic_block_partitioner.bin", f1.bins_ptr(),
                                          3);

  partitioner.append_pixels(pixels1.begin(), pixels1.end(), 50'000);
  partitioner.append_pixels(pixels2.begin(), pixels2.end(), 50'000);
  partitioner.finalize();

  std::size_t num_interactions = 0;
  for (const auto& [bid, _] : partitioner.block_index()) {
    const auto blk = partitioner.merge_blocks(bid);
    num_interactions += static_cast<std::size_t>(blk.nRecords);
  }

  CHECK(num_interactions == pixels1.size() + pixels2.size());
}

TEST_CASE("devel") {
  spdlog::default_logger()->set_level(spdlog::level::debug);
  // const std::vector<std::uint32_t> resolutions =
  //    hic::utils::list_resolutions((datadir / "4DNFIZ1ZVXC8.hic9").string());
  const std::vector<std::uint32_t> resolutions{500'000, 2'500'000};
  {
    const hic::File f1((datadir / "4DNFIZ1ZVXC8.hic9").string(), resolutions.front());
    // clang-format off
    const HiCHeader header{
            "/tmp/test.hic",           // url
            9,                         // version
            -1,                        // masterIndexOffset
            "hg38",                    // genomeID
            -1,                        // nviPosition
            -1,                        // nviLength
            f1.chromosomes(),
            {resolutions},         // resolutions
            {{"software", "hictk"}}    // attributes
    };
    // clang-format on

    HiCFileWriter w(header);
    w.write_header();

    const auto sel = f1.fetch();
    w.add_pixels(sel.begin<float>(), sel.end<float>());
    w.serialize();
  }

  // for (const auto& resolution : resolutions) {
  //   const hic::File f1((datadir / "4DNFIZ1ZVXC8.hic9").string(), resolution);
  //   const hic::File f2("/tmp/test.hic", resolution);
  //   const auto expected_pixels = f1.fetch().read_all<float>();
  //   const auto pixels = f2.fetch().read_all<float>();
  //
  //   REQUIRE(expected_pixels.size() == pixels.size());
  //   for (std::size_t i = 0; i < pixels.size(); ++i) {
  //     CHECK(expected_pixels[i] == pixels[i]);
  //   }
  // }
}

TEST_CASE("devel1") {
  spdlog::default_logger()->set_level(spdlog::level::debug);
  HiCFileReader f1("/tmp/4DNFIZ1ZVXC8.hic");
  auto footer1 = f1.read_footer(1, 1, MatrixType::observed, balancing::Method::NONE(),
                                MatrixUnit::BP, 2'500'000);
  auto idx1 = f1.read_index(footer1.fileOffset(), Chromosome{1, "chr2L", 1},
                            Chromosome{1, "chr2L", 1}, MatrixUnit::BP, 2'500'000);

  HiCFileReader f2("/tmp/test.hic");
  auto footer2 = f2.read_footer(1, 1, MatrixType::observed, balancing::Method::NONE(),
                                MatrixUnit::BP, 2'500'000);
  auto idx2 = f2.read_index(footer2.fileOffset(), Chromosome{1, "chr2L", 1},
                            Chromosome{1, "chr2L", 1}, MatrixUnit::BP, 2'500'000);
}

}  // namespace hictk::hic::test::file_writer
