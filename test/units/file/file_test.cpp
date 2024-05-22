// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/file.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <string>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/pixel.hpp"

using namespace hictk;

namespace hictk::test::file {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("File", "[file][short]") {
  const std::uint32_t resolution = 1'000'000;
  const auto path_hic = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();
  const auto path_cooler = (datadir / "integration_tests" / "4DNFIZ1ZVXC8.mcool").string();

  const auto uri_cooler = fmt::format(FMT_STRING("{}::/resolutions/{}"), path_cooler, resolution);

  const cooler::File ref(uri_cooler);

  SECTION("ctors") {
    CHECK(File(path_hic, resolution).path() == path_hic);
    CHECK(File(path_cooler, resolution).path() == path_cooler);
    CHECK(File(uri_cooler).uri() == uri_cooler);

    CHECK_THROWS(File(path_cooler, resolution, hic::MatrixType::expected));
    CHECK_THROWS(File(path_cooler, resolution, hic::MatrixType::observed, hic::MatrixUnit::FRAG));
  }

  SECTION("accessors") {
    SECTION("hic") {
      const hic::File ref_hic(path_hic, resolution);
      const auto hf = File(path_hic, resolution);

      CHECK(hf.is_hic());
      CHECK(hf.path() == path_hic);
      CHECK(hf.uri() == path_hic);

      CHECK(hf.chromosomes() == ref_hic.chromosomes());
      CHECK(hf.bins() == ref_hic.bins());

      CHECK(hf.resolution() == ref_hic.resolution());
      CHECK(hf.nbins() == ref_hic.nbins());
      CHECK(hf.nchroms() == ref_hic.nchroms());
    }

    SECTION("cooler") {
      const auto clr = File(path_cooler, resolution);

      CHECK(clr.is_cooler());
      CHECK(clr.path() == path_cooler);
      CHECK(clr.uri() == uri_cooler);

      CHECK(clr.chromosomes() == ref.chromosomes());
      CHECK(clr.bins() == ref.bins());

      CHECK(clr.resolution() == ref.resolution());
      CHECK(clr.nbins() == ref.nbins());
      CHECK(clr.nchroms() == ref.nchroms());
    }
  }

  SECTION("fetch") {
    SECTION("hic") {
      const auto hf = File(path_hic, resolution);
      const auto sel1 = ref.fetch("chr4", 0, 1'000'000);
      const auto sel2 = hf.fetch("chr4", 0, 1'000'000);

      auto first1 = sel1.begin<std::int32_t>();
      auto last1 = sel1.end<std::int32_t>();

      auto first2 = sel2.begin<std::int32_t>();
      auto last2 = sel2.end<std::int32_t>();
      CHECK(std::distance(first1, last1) == std::distance(first2, last2));
    }

    SECTION("hic gw") {
      const auto hf = File(path_hic, resolution);
      const auto sel1 = ref.fetch();
      const auto sel2 = hf.fetch();

      auto first1 = sel1.begin<std::int32_t>();
      auto last1 = sel1.end<std::int32_t>();

      auto first2 = sel2.begin<std::int32_t>();
      auto last2 = sel2.end<std::int32_t>();
      CHECK(std::distance(first1, last1) == std::distance(first2, last2));
    }

    SECTION("cooler") {
      const auto clr = File(path_cooler, resolution);
      const auto sel1 = ref.fetch("chr4");
      const auto sel2 = clr.fetch("chr4");

      auto first1 = sel1.begin<std::int32_t>();
      auto last1 = sel1.end<std::int32_t>();

      auto first2 = sel2.begin<std::int32_t>();
      auto last2 = sel2.end<std::int32_t>();
      CHECK(std::distance(first1, last1) == std::distance(first2, last2));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("PixelSelector", "[file][short]") {
  const std::uint32_t resolution = 1'000'000;
  const auto path_hic = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();
  const auto path_cooler = (datadir / "integration_tests" / "4DNFIZ1ZVXC8.mcool").string();

  SECTION("accessors") {
    SECTION("hic") {
      const auto hf = File(path_hic, resolution);
      auto sel1 = hf.fetch("chr2L", "chr2R");

      CHECK(sel1.coord1().bin1.chrom().name() == "chr2L");
      CHECK(sel1.coord2().bin1.chrom().name() == "chr2R");
      CHECK(sel1.bins().resolution() == resolution);

      CHECK(sel1.read_all<std::int32_t>().size() == 624);
    }
    SECTION("hic gw") {
      const auto hf = File(path_hic, resolution);
      auto sel1 = hf.fetch();
      CHECK(sel1.coord1() == PixelCoordinates{});
      CHECK(sel1.coord2() == PixelCoordinates{});

      CHECK(sel1.read_all<std::int32_t>().size() == 10'148);
    }
    SECTION("cooler") {
      const auto clr = File(path_cooler, resolution);
      auto sel1 = clr.fetch("chr2L", "chr2R");

      CHECK(sel1.coord1().bin1.chrom().name() == "chr2L");
      CHECK(sel1.coord2().bin1.chrom().name() == "chr2R");
      CHECK(sel1.bins().resolution() == resolution);

      CHECK(sel1.read_all<std::int32_t>().size() == 624);
    }
  }
}

}  // namespace hictk::test::file
