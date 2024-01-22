// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>
#include <vector>

#include "./tmpdir.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk::hic::test::utils {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiC: utils merge", "[merge][utils][long]") {
  SECTION("merge gw") {
    const auto src = datadir / "4DNFIZ1ZVXC8.hic9";
    const auto dest = testdir() / "hic_merge_test_001.hic";

    const std::uint32_t resolution = 500'000;
    const std::array<std::string, 2> sources{src.string(), src.string()};
    hic::utils::merge(sources.begin(), sources.end(), dest.string(), resolution, testdir(), true,
                      1'000);

    const File f1(src.string(), resolution);
    const File f2(dest.string(), resolution);

    const auto pixels1 = f1.fetch().read_all<float>();
    const auto pixels2 = f2.fetch().read_all<float>();

    REQUIRE(pixels1.size() == pixels2.size());
    for (std::size_t i = 0; i < pixels1.size(); ++i) {
      CHECK(pixels1[i].coords == pixels2[i].coords);
      CHECK(pixels1[i].count * 2 == pixels2[i].count);
    }
  }  // namespace hictk::hic::test::utils

  SECTION("merge chromosomes") {
    const auto src = datadir / "4DNFIZ1ZVXC8.hic9";
    const auto dest = testdir() / "hic_merge_test_002.hic";
    std::vector<std::string> sources{};
    const std::uint32_t resolution = 500'000;
    {
      spdlog::default_logger()->set_level(spdlog::level::warn);
      const File f(src.string(), resolution);

      for (std::uint32_t chrom1_id = 0; chrom1_id < f.chromosomes().size(); ++chrom1_id) {
        const auto& chrom1 = f.chromosomes().at(chrom1_id);
        if (chrom1.is_all()) {
          continue;
        }
        for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < f.chromosomes().size(); ++chrom2_id) {
          const auto& chrom2 = f.chromosomes().at(chrom2_id);

          const auto sel = f.fetch(chrom1.name(), chrom2.name());
          if (sel.empty()) {
            continue;
          }

          sources.emplace_back((testdir() / fmt::format(FMT_STRING("hic_merge_test_002.{}_{}.hic"),
                                                        chrom1.name(), chrom2.name()))
                                   .string());

          hic::internal::HiCFileWriter w(sources.back(), f.chromosomes(), {f.bin_size()}, "", 1,
                                         1'000, testdir());

          w.add_pixels(resolution, sel.begin<float>(), sel.end<float>());
          w.serialize();
        }
      }
    }

    spdlog::default_logger()->set_level(spdlog::level::info);
    hic::utils::merge(sources.begin(), sources.end(), dest.string(), resolution, testdir(), true,
                      1'000);

    const File f1(src.string(), resolution);
    const File f2(dest.string(), resolution);

    const auto pixels1 = f1.fetch().read_all<float>();
    const auto pixels2 = f2.fetch().read_all<float>();

    REQUIRE(pixels1.size() == pixels2.size());
    for (std::size_t i = 0; i < pixels1.size(); ++i) {
      CHECK(pixels1[i] == pixels2[i]);
    }
  }

  SECTION("merge - different reference") {
    const auto src1 = datadir / "4DNFIZ1ZVXC8.hic9";
    const auto src2 = datadir / "ENCFF993FGR.2500000.hic";
    const auto dest = testdir() / "cooler_merge_test_003.cool";

    const std::array<std::string, 2> sources{src1.string(), src2.string()};

    CHECK_THROWS_WITH(hic::utils::merge(sources.begin(), sources.end(), dest.string(), 2'500'000),
                      Catch::Matchers::ContainsSubstring("use different reference genomes"));
  }
}

}  // namespace hictk::hic::test::utils
