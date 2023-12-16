// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <catch2/catch_test_macros.hpp>
#include <filesystem>

#include "hictk/cooler/validation.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: format checking", "[cooler][short]") {
  SECTION("test .cool") {
    const auto path = datadir / "cooler_test_file.cool";
    CHECK(utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
  }

  SECTION("test .mcool") {
    const auto path = datadir / "multires_cooler_test_file.mcool";
    constexpr auto suffix{"::/resolutions/400000"};

    CHECK(!utils::is_cooler(path.string()));
    CHECK(utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
    CHECK(utils::is_cooler(path.string() + suffix));
  }

  SECTION("test .scool") {
    const auto path = datadir / "single_cell_cooler_test_file.scool";
    constexpr auto suffix{"::/cells/GSM2687248_41669_ACAGTG-R1-DpnII.100000.cool"};

    CHECK(!utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(utils::is_scool_file(path.string()));
    CHECK(utils::is_cooler(path.string() + suffix));
  }

  SECTION("test with empty .h5 file") {
    const auto path = datadir / "empty_test_file.h5";
    CHECK(!utils::is_cooler(path.string()));
    CHECK(!utils::is_multires_file(path.string()));
    CHECK(!utils::is_scool_file(path.string()));
  }

  SECTION("test with nonexistent file") {
    const auto invalid_path = datadir / "void.nonexistent";
    CHECK(utils::is_cooler(invalid_path.string()).unable_to_open_file);
    CHECK(utils::is_multires_file(invalid_path.string()).unable_to_open_file);
    CHECK(utils::is_scool_file(invalid_path.string()).unable_to_open_file);
  }

  SECTION("test corrupted .cool") {
    SECTION("missing format attribute") {
      const auto path = datadir / "invalid_coolers/missing_format_attr.cool";
      CHECK(utils::is_cooler(path.string()).missing_or_invalid_format_attr);
    }
    SECTION("invalid format attribute") {
      const auto path = datadir / "invalid_coolers/invalid_format_attr.cool";
      CHECK(utils::is_cooler(path.string()).missing_or_invalid_format_attr);
    }
  }

  SECTION("test corrupted .mcool") {
    // This file is missing group /resolutions/400000/pixels
    const auto path = datadir / "invalid_coolers/missing_pixels_group.mcool";
    const auto status = utils::is_multires_file(path.string());

    CHECK(!status);
    CHECK(status.is_hdf5);
    CHECK(!status.is_multires_file);
    CHECK(!status.missing_or_invalid_format_attr);
    CHECK(!status.missing_or_invalid_bin_type_attr);
    CHECK(status.uri == path.string());
    CHECK(status.missing_groups.empty());

    REQUIRE(status.invalid_resolutions.size() == 1);
    const auto& invalid_res = status.invalid_resolutions.front();

    const auto corrupted_uri_expected =
        fmt::format(FMT_STRING("{}::/resolutions/400000"), path.string());
    CHECK(invalid_res.uri == corrupted_uri_expected);
    CHECK(!invalid_res.is_cooler);
    REQUIRE(invalid_res.missing_groups.size() == 1);
    CHECK(invalid_res.missing_groups.front() == "pixels");
  }

  SECTION("test corrupted .scool") {
    // In this file, the number of groups under /cells and number of cells from ncells attribute
    // mismatch
    const auto path = datadir / "invalid_coolers/invalid_ncells_attribute.scool";
    const auto status = utils::is_scool_file(path.string());

    CHECK(!status);
    CHECK(status.is_hdf5);
    CHECK(!status.is_scool_file);
    CHECK(!status.missing_or_invalid_format_attr);
    CHECK(!status.missing_or_invalid_bin_type_attr);
    CHECK(status.uri == path.string());
    CHECK(status.missing_groups.empty());
    CHECK(status.unexpected_number_of_cells);
    CHECK(status.invalid_cells.empty());
  }
}

TEST_CASE("Cooler: index validation", "[cooler][short]") {
  SECTION("valid index") {
    const auto path1 = datadir / "ENCFF993FGR.2500000.cool";
    const auto path2 = datadir / "cooler_test_file.cool";
    CHECK(cooler::utils::index_is_valid(path1.string()));
    CHECK(cooler::utils::index_is_valid(path2.string()));
  }
  SECTION("broken index") {
    SECTION(".cool") {
      const auto path = datadir / "invalid_coolers/corrupted_index.mcool::/resolutions/10000000";
      CHECK_FALSE(cooler::utils::index_is_valid(path.string()));
    }
    SECTION(".mcool") {
      const auto path = datadir / "invalid_coolers/corrupted_index.mcool";
      CHECK_FALSE(cooler::utils::index_is_valid(path.string()));
    }
  }
}

}  // namespace hictk::cooler::test::cooler_file
