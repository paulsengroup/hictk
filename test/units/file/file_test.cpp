// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/file.hpp"

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <string>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/pixel.hpp"
#include "hictk/string.hpp"
#include "hictk/test/testdir.hpp"

using namespace hictk;

namespace hictk::test::file {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("File", "[file][short]") {
  const std::uint32_t resolution = 1'000'000;
  const auto path_hic = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();
  const auto path_cooler = (datadir / "cooler" / "4DNFIZ1ZVXC8.mcool").string();

  const auto uri_cooler = fmt::format(FMT_STRING("{}::/resolutions/{}"), path_cooler, resolution);

  const cooler::File ref(uri_cooler);

  SECTION("ctors") {
    const auto path_singleres_hic = (datadir / "hic" / "ENCFF993FGR.2500000.hic").string();
    const auto path_singleres_mcool =
        (datadir / "cooler" / "singleres_cooler_test_file.mcool").string();
    SECTION("valid") {
      CHECK(File(path_hic, resolution).path() == path_hic);
      CHECK(File(path_cooler, resolution).path() == path_cooler);
      CHECK(File(uri_cooler).uri() == uri_cooler);
      CHECK(File(path_singleres_hic).resolution() == 2'500'000);
      CHECK(File(path_singleres_mcool).resolution() == 6'400'000);
    }

    SECTION("invalid") {
      // Invalid params for .cool file
      CHECK_THROWS_WITH(
          File(fmt::format(FMT_STRING("{}::/resolutions/{}"), path_cooler, resolution + 1)),
          Catch::Matchers::ContainsSubstring("resolution is required"));
      CHECK_THROWS_WITH(File(path_cooler, resolution + 1),
                        Catch::Matchers::ContainsSubstring("unable to find resolution"));
      CHECK_THROWS_WITH(File(uri_cooler, resolution + 1),
                        Catch::Matchers::ContainsSubstring("found an unexpected resolution"));

      // Invalid params for .mcool files
      CHECK_THROWS_WITH(File(path_cooler),
                        Catch::Matchers::ContainsSubstring("resolution is required"));
      CHECK_THROWS_WITH(File(path_cooler, resolution, hic::MatrixType::expected),
                        Catch::Matchers::ContainsSubstring("should always be \"observed\""));
      CHECK_THROWS_WITH(
          File(path_cooler, resolution, hic::MatrixType::observed, hic::MatrixUnit::FRAG),
          Catch::Matchers::ContainsSubstring("should always be \"BP\""));

      // Invalid params for .hic files
      CHECK_THROWS_WITH(File(path_hic),
                        Catch::Matchers::ContainsSubstring("resolution is required"));
    }
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
      CHECK(hf.nchroms(true) == ref_hic.nchroms(true));
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
      CHECK(clr.nchroms(true) == ref.nchroms());
    }
  }

  SECTION("fetch") {
    SECTION("hic") {
      const auto hf = File(path_hic, resolution);
      const auto sel1 = ref.fetch("chr4", 0, 1'000'000);
      const auto sel2 = hf.fetch("chr4", 0, 1'000'000);
      CHECK(sel1.size() == sel2.size());

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
      CHECK(sel1.size() == sel2.size());

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
      CHECK(sel1.size() == sel2.size());

      auto first1 = sel1.begin<std::int32_t>();
      auto last1 = sel1.end<std::int32_t>();

      auto first2 = sel2.begin<std::int32_t>();
      auto last2 = sel2.end<std::int32_t>();
      CHECK(std::distance(first1, last1) == std::distance(first2, last2));
    }
  }
}

TEST_CASE("PixelSelector", "[file][short]") {
  const std::uint32_t resolution = 1'000'000;
  const auto path_hic = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();
  const auto path_cooler = (datadir / "cooler" / "4DNFIZ1ZVXC8.mcool").string();

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

TEST_CASE("PixelSelector::iterator", "[file][short]") {
  const std::uint32_t resolution = 1'000'000;
  const auto hf =
      std::make_shared<const File>((datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string(), resolution);
  const auto clr = std::make_shared<const File>(
      (datadir / "cooler" / "4DNFIZ1ZVXC8.mcool").string(), resolution);

  const std::array<
      std::tuple<std::string, std::shared_ptr<const File>, std::shared_ptr<const File>>, 3>
      files{
          // clang-format off
          std::make_tuple("hic", hf, hf),
          std::make_tuple("hic gw", hf, hf),
          std::make_tuple("cooler", clr, clr),
          // clang-format on
      };

  using T = std::int32_t;

  for (const auto& [label, f1, f2] : files) {
    SECTION(label) {
      const auto is_gw = internal::ends_with(label, "gw");
      const auto sel1 = is_gw ? f1->fetch() : f1->fetch("chr2L", "chr2R");
      const auto sel2 = is_gw ? f2->fetch() : f2->fetch("chr2L", "chr2R");

      SECTION("operator==") {
        auto it1 = sel1.begin<T>();
        auto it2 = sel2.begin<T>();
        CHECK(it1 == it2);

        it1 = ++sel1.begin<T>();
        it2 = ++sel2.begin<T>();
        CHECK(it1 == it2);

        it1 = sel1.end<T>();
        it2 = sel2.end<T>();
        CHECK(it1 == it2);
      }

      SECTION("operator!=") {
        auto it1 = sel1.begin<T>();
        auto it2 = ++sel2.begin<T>();
        CHECK(it1 != it2);

        it2 = sel2.end<T>();
        CHECK(it1 != it2);
      }
    }
  }
  SECTION("cooler-hic") {
    const auto sel1 = clr->fetch("chr2L", "chr2R");
    const auto sel2 = hf->fetch("chr2L", "chr2R");

    SECTION("operator==") {
      auto it1 = sel1.begin<T>();
      auto it2 = sel2.begin<T>();
      CHECK_FALSE(it1 == it2);

      it1 = ++sel1.begin<T>();
      it2 = ++sel2.begin<T>();
      CHECK_FALSE(it1 == it2);

      it1 = sel1.end<T>();
      it2 = sel2.end<T>();
      CHECK_FALSE(it1 == it2);
    }

    SECTION("operator!=") {
      auto it1 = sel1.begin<T>();
      auto it2 = ++sel2.begin<T>();
      CHECK(it1 != it2);

      it2 = sel2.end<T>();
      CHECK(it1 != it2);
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::file
