// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>

#include "hictk/balancing/methods.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/test/testdir.hpp"

using namespace hictk::hic;

namespace hictk::hic::test::file {

static const auto& datadir = hictk::test::datadir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)

// NOLINTBEGIN(cert-err58-cpp)
const auto single_res = (datadir / "hic" / "ENCFF993FGR.2500000.hic").string();
const auto pathV8 = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();
const auto pathV9 = (datadir / "hic" / "4DNFIZ1ZVXC8.hic9").string();
const auto path_binary = (datadir / "various" / "data.zip").string();
// NOLINTEND(cert-err58-cpp)

TEST_CASE("HiC: utils is_hic_file", "[hic][short]") {
  CHECK(utils::is_hic_file(pathV8));
  CHECK_FALSE(utils::is_hic_file(path_binary));
}

TEST_CASE("HiC: ctors", "[hic][short]") {
  SECTION("valid") {
    CHECK_NOTHROW(File{pathV8, 1'000});
    CHECK_NOTHROW(File{pathV9, 1'000});
    CHECK(File{single_res}.resolution() == 2'500'000);
  }
  SECTION("invalid") {
    CHECK_THROWS_WITH(File(pathV8, 1), Catch::Matchers::ContainsSubstring(
                                           "does not have interactions for resolution"));
    CHECK_THROWS_WITH(File(pathV8), Catch::Matchers::ContainsSubstring("resolution is required"));
  }
}

TEST_CASE("HiC: file accessors", "[hic][short]") {
  File f(pathV8, 1'000);

  CHECK(f.path() == pathV8);
  CHECK(f.name() == pathV8);
  CHECK(f.version() == 8);
  CHECK(f.chromosomes().size() == 9);
  CHECK(f.nchroms(true) == 9);
  CHECK(f.nchroms() == 8);
  CHECK(f.nbins() == 137572);
  CHECK(f.assembly() == "dm6");
  CHECK(f.attributes().at("software") == "Juicer Tools Version 1.22.01");

  CHECK(f.avail_resolutions().size() == 10);
  CHECK(f.avail_resolutions().front() == 1'000);
  CHECK(f.avail_resolutions().back() == 2'500'000);

  CHECK(f.avail_normalizations().size() == 4);
  CHECK(f.avail_normalizations().front() == "KR");
  CHECK(f.avail_normalizations().back() == "VC_SQRT");

  CHECK(f.open(2'500'000).resolution() == 2'500'000);

  SECTION("invalid") {
    CHECK_THROWS(File(pathV8, std::numeric_limits<std::uint32_t>::max(), MatrixType::observed,
                      MatrixUnit::BP));
    CHECK_THROWS(File("non-existing-file", 1));
    CHECK_THROWS(File("https://localhost:non-existing-url", 1));
    CHECK_THROWS(File("test/CMakeLists.txt", 1));
  }
}

TEST_CASE("HiC: footer cache", "[hic][short]") {
  File f(pathV8, 2'500'000, MatrixType::observed, MatrixUnit::BP, 1);

  CHECK(f.num_cached_footers() == 0);
  for (const auto& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    std::ignore = f.fetch(chrom.name());
  }

  CHECK(f.num_cached_footers() == 8);

  const auto sel1 = f.fetch("chr2L");
  const auto sel2 = f.fetch("chr2L");

  // this check relies on the fact that metadata are stored in footers, and that footers are
  // looked up in the cache when creating matrix selectors
  CHECK(&sel1.metadata() == &sel2.metadata());

  f.purge_footer_cache();
  CHECK(f.num_cached_footers() == 0);

  const auto sel3 = f.fetch("chr2L");

  CHECK(f.num_cached_footers() == 1);
  CHECK(&sel1.metadata() != &sel3.metadata());
}

TEST_CASE("HiC: fetch", "[hic][short]") {
  const auto norm = hictk::balancing::Method::NONE();
  const File f(pathV8, 2'500'000, MatrixType::observed, MatrixUnit::BP);

  REQUIRE(f.chromosomes().size() == 9);

  const auto* chrom1 = "chr2L";
  const auto* chrom2 = "chr2R";
  SECTION("intra-chromosomal") {
    SECTION("valid") {
      auto sel = f.fetch(chrom1, norm);
      CHECK(sel.chrom1() == chrom1);
    }
    SECTION("invalid") {
      CHECK_THROWS_WITH(f.fetch("chr2L:2,500,000-5,000,000", "chr2L:0-2,500,000"),
                        Catch::Matchers::ContainsSubstring("overlaps with the lower-triangle"));
    }
  }

  SECTION("inter-chromosomal") {
    SECTION("valid") {
      auto sel = f.fetch(chrom1, chrom2, norm);
      CHECK(sel.chrom1() == chrom1);
      CHECK(sel.chrom2() == chrom2);
    }
    SECTION("invalid") {
      SECTION("invalid") {
        CHECK_THROWS_WITH(f.fetch(chrom2, chrom1),
                          Catch::Matchers::ContainsSubstring("overlaps with the lower-triangle"));
      }
    }
  }

  SECTION("valid, but empty matrix") {
    auto sel = f.fetch("chrM", norm);
    const auto buff = sel.read_all<float>();
    CHECK(buff.empty());
  }

  SECTION("invalid chromosome") {
    CHECK_THROWS(f.fetch("not-a-chromosome", norm));
    CHECK_THROWS(f.fetch("chr2L", "not-a-chromosome", norm));
  }

  SECTION("malformed") {
    CHECK_THROWS(f.fetch(chrom2, chrom1, norm));  // NOLINT

    // Matrix does not have contacts for fragments
    CHECK_THROWS(
        File(pathV8, f.resolution(), MatrixType::observed, MatrixUnit::FRAG).fetch(chrom1, norm));
  }
}

TEST_CASE("HiC: list_resolutions", "[hic][short]") {
  constexpr std::array<std::uint32_t, 10> expected{1000,   5000,   10000,  25000,   50000,
                                                   100000, 250000, 500000, 1000000, 2500000};

  auto found = utils::list_resolutions(pathV8);
  REQUIRE(found.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    CHECK(expected[i] == found[i]);
  }

  found = utils::list_resolutions(pathV9);
  REQUIRE(found.size() == expected.size());
  for (std::size_t i = 0; i < expected.size(); ++i) {
    CHECK(expected[i] == found[i]);
  }
}

TEST_CASE("HiC: list_normalizations", "[hic][short]") {
  const std::array<balancing::Method, 4> expected_v8{
      balancing::Method::KR(),
      balancing::Method::SCALE(),
      balancing::Method::VC(),
      balancing::Method::VC_SQRT(),
  };

  const std::array<balancing::Method, 3> expected_v9{
      balancing::Method::SCALE(),
      balancing::Method::VC(),
      balancing::Method::VC_SQRT(),
  };

  for (const auto* policy : {"union", "intersection"}) {
    auto found = utils::list_normalizations(pathV8, policy);
    REQUIRE(found.size() == expected_v8.size());
    for (std::size_t i = 0; i < expected_v8.size(); ++i) {
      CHECK(expected_v8[i] == found[i]);
    }

    found = utils::list_normalizations(pathV9, policy);
    REQUIRE(found.size() == expected_v9.size());
    for (std::size_t i = 0; i < expected_v9.size(); ++i) {
      CHECK(expected_v9[i] == found[i]);
    }
  }

  CHECK_THROWS_AS(utils::list_normalizations(pathV8, "invalid"), std::invalid_argument);
}
// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::hic::test::file
