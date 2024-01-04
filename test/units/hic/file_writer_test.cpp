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
#include "hictk/hic/file_reader.hpp"
#include "hictk/reference.hpp"
#include "tmpdir.hpp"

using namespace hictk::hic;

namespace hictk::hic::test::file_writer {

using namespace hictk::hic::internal;

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
    w.write_master_index_offset(master_index_offset);

    std::ofstream ofs(header1.url, std::ios::binary | std::ios::app);
    ofs.write("0", 1);  // Add padding
  }

  const auto header2 = HiCFileReader{header1.url}.header();
  CHECK(header1 == header2);
}

TEST_CASE("HiC: BlockMapperIntra", "[hic][v9][short]") {
  SECTION("intra") {
    // Test case based on blocks fetched by running an instrumented version of
    // hictk dump test/data/hic/4DNFIZ1ZVXC8.hic9 --resolution 10000 --range chr3R:0-50000
    HiCFileWriter::BlockMapperIntra mapper(803, 4);
    CHECK(mapper(0, 0) == 0);
    CHECK(mapper(0, 100) == 0);
    CHECK(mapper(802, 802) == 0);
    CHECK(mapper(803, 803) == 1);
    CHECK(mapper(1038, 2137) == 1);
    CHECK(mapper(235, 1376) == 5);
    CHECK(mapper(8, 3203) == 5);
  }

  SECTION("inter") {
    // Test case based on blocks fetched by running an instrumented version of
    // hictk dump test/data/hic/4DNFIZ1ZVXC8.hic9 --resolution 10000 --range chr3R:0-50000 --range2
    // chr3R:0-10000000
    HiCFileWriter::BlockMapperInter mapper(803, 4);
    CHECK(mapper(0, 0) == 0);
    CHECK(mapper(0, 100) == 0);
    CHECK(mapper(802, 802) == 0);
    CHECK(mapper(7, 803) == 4);
    CHECK(mapper(795, 1605) == 4);
  }
}

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

    const std::vector<ThinPixel<float>> pixels{{0, 0, 1}, {0, 1, 2}, {1, 1, 3}};
    const InteractionBlock blk(0, 2, pixels);

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
        w.write_interaction_block(blk, 0, 0);
      }

      const auto [master_index_offset, matrix_metadata_size] =
          w.write_body_metadata(chrom.id(), chrom.id());
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

}  // namespace hictk::hic::test::file_writer
