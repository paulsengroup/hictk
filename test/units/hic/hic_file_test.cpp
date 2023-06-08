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
  HiCFile f(pathV8);

  CHECK(f.url() == pathV8);
  CHECK(f.name() == pathV8);
  CHECK(f.version() == 8);
  CHECK(f.chromosomes().size() == 9);
  CHECK(f.assembly() == "dm6");

  CHECK(f.resolutions().size() == 10);
  CHECK(f.resolutions().front() == 2500000);
  CHECK(f.resolutions().back() == 1000);

  SECTION("invalid") {
    CHECK_THROWS(HiCFile("non-existing-file"));
    CHECK_THROWS(HiCFile("https://localhost:non-existing-url"));
    CHECK_THROWS(HiCFile("test/CMakeLists.txt"));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile footer cache", "[hic][short]") {
  HiCFile f(pathV8);

  REQUIRE(f.resolutions().size() == 10);

  CHECK(f.num_cached_footers() == 0);
  for (const auto res : f.resolutions()) {
    std::ignore = f.get_matrix_selector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                        MatrixUnit::BP, res);
  }

  CHECK(f.num_cached_footers() == f.resolutions().size());

  const auto sel1 = f.get_matrix_selector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                          MatrixUnit::BP, 2500000);
  const auto sel2 = f.get_matrix_selector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                          MatrixUnit::BP, 2500000);

  // this check relies on the fact that chrom1Norm are stored in the footer, and that footers are
  // looked up in the cache when creating matrix selectors
  CHECK(&sel1.chrom1Norm() == &sel2.chrom1Norm());

  f.purge_footer_cache();
  CHECK(f.num_cached_footers() == 0);

  const auto sel3 = f.get_matrix_selector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                          MatrixUnit::BP, 2500000);

  CHECK(f.num_cached_footers() == 1);
  CHECK(&sel1.chrom1Norm() != &sel3.chrom1Norm());
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile get_matrix_selector", "[hic][short]") {
  HiCFile f(pathV8);

  REQUIRE(f.chromosomes().size() == 9);

  constexpr auto mt = MatrixType::observed;
  constexpr auto norm = NormalizationMethod::NONE;
  constexpr auto unit = MatrixUnit::BP;
  constexpr std::int32_t res = 2500000;

  const auto chrom1 = f.chromosomes().at("chr2L");
  const auto chrom2 = f.chromosomes().at("chr2R");

  SECTION("intra-chromosomal") {
    auto sel = f.get_matrix_selector(chrom1.name, mt, norm, unit, res);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.isIntra());

    sel = f.get_matrix_selector(chrom1.index, mt, norm, unit, res);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.isIntra());

    sel = f.get_matrix_selector(chrom1, mt, norm, unit, res);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.isIntra());
  }

  SECTION("inter-chromosomal") {
    auto sel = f.get_matrix_selector(chrom1.name, chrom2.name, mt, norm, unit, res);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.chrom2() == chrom2);

    sel = f.get_matrix_selector(chrom1.index, chrom2.index, mt, norm, unit, res);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.chrom2() == chrom2);

    sel = f.get_matrix_selector(chrom1, chrom2, mt, norm, unit, res);
    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.chrom2() == chrom2);

    CHECK(sel.chrom1() == chrom1);
    CHECK(sel.chrom2() == chrom2);
  }

  SECTION("valid, but empty matrix") {
    auto sel = f.get_matrix_selector("chrM", mt, norm, unit, res);
    std::vector<contactRecord> buff{};
    sel.fetch(buff);
    CHECK(buff.empty());
  }

  SECTION("invalid chromosome") {
    CHECK_THROWS(f.get_matrix_selector("not-a-chromosome", mt, norm, unit, res));
    CHECK_THROWS(f.get_matrix_selector(chrom1.name, "not-a-chromosome", mt, norm, unit, res));
    CHECK_THROWS(f.get_matrix_selector(999, mt, norm, unit, res));
    CHECK_THROWS(f.get_matrix_selector(chrom1.index, 999, mt, norm, unit, res));
  }

  SECTION("malformed") {
    CHECK_THROWS(f.get_matrix_selector(chrom2, chrom1, mt, norm, unit, res));  // NOLINT
    CHECK_THROWS(f.get_matrix_selector(chrom1, mt, norm, unit, 123));
    CHECK_THROWS(
        f.get_matrix_selector(chrom1, MatrixType::expected, NormalizationMethod::VC, unit, res));

    // Matrix does not have contacts for fragments
    CHECK_THROWS(f.get_matrix_selector(chrom1, mt, norm, MatrixUnit::FRAG, res));
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("GenomicCoordinates", "[hic][short]") {
  SECTION("chrom-only") {
    const auto coord = GenomicCoordinates::fromString("chr1");
    CHECK(coord.chrom == "chr1");
    CHECK(coord.start == 0);
    CHECK(coord.end == 0);
  }

  SECTION("valid") {
    auto coord = GenomicCoordinates::fromString("chr1:0-1000");
    CHECK(coord.chrom == "chr1");
    CHECK(coord.start == 0);
    CHECK(coord.end == 1000);

    coord = GenomicCoordinates::fromString("chr1:0:1000");
    CHECK(coord.chrom == "chr1");
    CHECK(coord.start == 0);
    CHECK(coord.end == 1000);
  }

  SECTION("invalid") {
    CHECK_THROWS(GenomicCoordinates::fromString("chr1:0"));
    CHECK_THROWS(GenomicCoordinates::fromString("chr1:0:"));

    CHECK_THROWS(GenomicCoordinates::fromString("chr1:0-"));
    CHECK_THROWS(GenomicCoordinates::fromString("chr1:-"));

    CHECK_THROWS(GenomicCoordinates::fromString("chr1::"));
    CHECK_THROWS(GenomicCoordinates::fromString("chr1:a:b"));
    CHECK_THROWS(GenomicCoordinates::fromString("chr1:a-b"));

    CHECK_THROWS(GenomicCoordinates::fromString("chr1:100-0"));
    CHECK_THROWS(GenomicCoordinates::fromString("chr1:0-100suffix"));
  }
}
