// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/multires_file.hpp"

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/test/testdir.hpp"

using namespace hictk;

namespace hictk::test::file {

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("MultiResFile", "[file][short]") {
  const auto path_hic = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();
  const auto path_mcool = (datadir / "cooler" / "4DNFIZ1ZVXC8.mcool").string();
  const std::uint32_t resolution = 1'000'000;

  SECTION("ctors") {
    CHECK(MultiResFile{path_hic}.path() == path_hic);
    CHECK(MultiResFile{hic::File(path_hic, resolution)}.path() == path_hic);

    CHECK(MultiResFile{path_mcool}.path() == path_mcool);
    CHECK(MultiResFile{cooler::MultiResFile(path_mcool, resolution)}.path() == path_mcool);

    CHECK_THROWS(MultiResFile{path_mcool, hic::MatrixType::expected});
    CHECK_THROWS(MultiResFile{path_mcool, hic::MatrixType::observed, hic::MatrixUnit::FRAG});
  }

  SECTION("accessors") {
    SECTION("hic") {
      CHECK(MultiResFile{path_hic}.is_hic());

      CHECK(MultiResFile{path_hic, hic::MatrixType::expected}.matrix_type() ==
            hic::MatrixType::expected);
      CHECK(MultiResFile{path_hic}.matrix_unit() == hic::MatrixUnit::BP);

      CHECK(MultiResFile{path_hic}.format() == "HIC");
      CHECK(MultiResFile{path_hic}.version() == 8);
      CHECK(MultiResFile{path_hic}.bin_type() == BinTable::Type::fixed);

      CHECK(MultiResFile{path_hic}.resolutions().size() == 10);
      CHECK(MultiResFile{path_hic}.chromosomes().size() == 9);
      CHECK(MultiResFile{path_hic}.avail_normalizations().size() == 4);
    }
    SECTION("mcool") {
      CHECK(MultiResFile{path_mcool}.is_mcool());

      CHECK(MultiResFile{path_mcool}.matrix_type() == hic::MatrixType::observed);
      CHECK(MultiResFile{path_mcool}.matrix_unit() == hic::MatrixUnit::BP);

      CHECK(MultiResFile{path_mcool}.format() == "HDF5::MCOOL");

      CHECK(MultiResFile{path_mcool}.version() == 2);
      CHECK(MultiResFile{path_mcool}.bin_type() == BinTable::Type::fixed);

      CHECK(MultiResFile{path_mcool}.resolutions().size() == 10);
      CHECK(MultiResFile{path_mcool}.chromosomes().size() == 8);
      CHECK(MultiResFile{path_mcool}.avail_normalizations().size() == 5);
    }
  }

  SECTION("open") {
    SECTION("hic") {
      CHECK(MultiResFile{path_hic}.open(resolution).resolution() == resolution);
      CHECK_THROWS(MultiResFile{path_hic}.open(resolution + 1));
    }

    SECTION("mcool") {
      CHECK(MultiResFile{path_mcool}.open(resolution).resolution() == resolution);
      CHECK_THROWS(MultiResFile{path_mcool}.open(resolution + 1));
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::file
