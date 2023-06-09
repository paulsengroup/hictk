// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <filesystem>
#include <string>

#include "hictk/hic.hpp"

using namespace hictk;

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

const auto pathV8 = (test::datadir / "4DNFIZ1ZVXC8.hic8").string();  // NOLINT(cert-err58-cpp)
const auto path_binary = (test::datadir / "data.zip").string();      // NOLINT(cert-err58-cpp)

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("utils: is_hic_file", "[hic][short]") {
  CHECK(utils::is_hic_file(pathV8));
  CHECK_FALSE(utils::is_hic_file(path_binary));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile accessors", "[hic][short]") {
  const HiCFile f(pathV8, 1'000);

  CHECK(f.url() == pathV8);
  CHECK(f.name() == pathV8);
  CHECK(f.version() == 8);
  CHECK(f.chromosomes().size() == 9);
  CHECK(f.assembly() == "dm6");

  CHECK(f.avail_resolutions().size() == 10);
  CHECK(f.avail_resolutions().front() == 2'500'000);
  CHECK(f.avail_resolutions().back() == 1000);

  CHECK(f.open_resolution(2'500'000).resolution() == 2'500'000);

  SECTION("invalid") {
    CHECK_THROWS(HiCFile(pathV8, (std::numeric_limits<std::uint32_t>::max)(), MatrixType::observed,
                         MatrixUnit::BP));
    CHECK_THROWS(HiCFile("non-existing-file", 1));
    CHECK_THROWS(HiCFile("https://localhost:non-existing-url", 1));
    CHECK_THROWS(HiCFile("test/CMakeLists.txt", 1));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile footer cache", "[hic][short]") {
  HiCFile f(pathV8, 2'500'000);

  CHECK(f.num_cached_footers() == 0);
  for (const auto& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    std::ignore = f.get_matrix_selector(chrom, NormalizationMethod::NONE);
  }

  CHECK(f.num_cached_footers() == 7);

  const auto sel1 = f.get_matrix_selector("chr2L", NormalizationMethod::NONE);
  const auto sel2 = f.get_matrix_selector("chr2L", NormalizationMethod::NONE);

  // this check relies on the fact that chrom1Norm are stored in the footer, and that footers are
  // looked up in the cache when creating matrix selectors
  CHECK(&sel1.chrom1Norm() == &sel2.chrom1Norm());

  f.purge_footer_cache();
  CHECK(f.num_cached_footers() == 0);

  const auto sel3 = f.get_matrix_selector("chr2L", NormalizationMethod::NONE);

  CHECK(f.num_cached_footers() == 1);
  CHECK(&sel1.chrom1Norm() != &sel3.chrom1Norm());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile get_matrix_selector", "[hic][short]") {
  constexpr auto norm = NormalizationMethod::NONE;
  HiCFile f(pathV8, 2'500'000, MatrixType::observed, MatrixUnit::BP);

  REQUIRE(f.chromosomes().size() == 9);

  const auto chrom1 = f.chromosomes().at("chr2L");
  const auto chrom2 = f.chromosomes().at("chr2R");
  SECTION("intra-chromosomal") {
    auto sel = f.get_matrix_selector(chrom1, norm);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.isIntra());

    sel = f.get_matrix_selector(chrom1.id(), norm);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.isIntra());

    sel = f.get_matrix_selector(chrom1, norm);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.isIntra());
  }

  SECTION("inter-chromosomal") {
    auto sel = f.get_matrix_selector(chrom1, chrom2, norm);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.chrom2() == chrom2);

    sel = f.get_matrix_selector(chrom1.id(), chrom2.id(), norm);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.chrom2() == chrom2);

    sel = f.get_matrix_selector(chrom1, chrom2, norm);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.chrom2() == chrom2);

    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.chrom2() == chrom2);
  }

  SECTION("valid, but empty matrix") {
    auto sel = f.get_matrix_selector("chrM", norm);
    std::vector<Pixel<float>> buff{};
    sel.fetch(buff);
    CHECK(buff.empty());
  }

  SECTION("invalid chromosome") {
    CHECK_THROWS(f.get_matrix_selector("not-a-chromosome", norm));
    CHECK_THROWS(f.get_matrix_selector(std::string{chrom1.name()}, "not-a-chromosome", norm));
    CHECK_THROWS(f.get_matrix_selector(999, norm));
    CHECK_THROWS(f.get_matrix_selector(chrom1.id(), 999, norm));
  }

  SECTION("malformed") {
    CHECK_THROWS(f.get_matrix_selector(chrom2, chrom1, norm));  // NOLINT
    CHECK_THROWS(HiCFile(pathV8, f.resolution(), MatrixType::expected, MatrixUnit::BP)
                     .get_matrix_selector(chrom1, NormalizationMethod::VC));

    // Matrix does not have contacts for fragments
    CHECK_THROWS(HiCFile(pathV8, f.resolution(), MatrixType::observed, MatrixUnit::FRAG)
                     .get_matrix_selector(chrom1, norm));
  }
}
