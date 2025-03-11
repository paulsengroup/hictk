// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/multires_cooler.hpp"

#include <fmt/format.h>

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <tuple>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/reference.hpp"
#include "hictk/test/testdir.hpp"

namespace hictk::cooler::test::multires_cooler_file {

static const auto& datadir = hictk::test::datadir;
static const auto& testdir = hictk::test::testdir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("MultiResCooler: open read-only", "[cooler][short]") {
  const auto path = datadir / "cooler" / "multires_cooler_test_file.mcool";

  auto mclr = MultiResFile(path.string());

  CHECK(mclr.chromosomes().size() == 20);
  CHECK(mclr.chromosomes().begin()->name() == "1");
  CHECK(mclr.resolutions().size() == 7);
  CHECK(mclr.attributes().format == MCOOL_MAGIC);
  CHECK(mclr.attributes().format_version == 2);
  CHECK(mclr.attributes().bin_type == BinTable::Type::fixed);

  CHECK(utils::is_cooler(mclr.open(1600000).uri()));
}

TEST_CASE("MultiResCooler: init files", "[cooler][short]") {
  const auto base_path = datadir / "cooler" / "cooler_test_file.cool";
  const auto base_resolution = File(base_path.string()).resolution();

  const auto path = testdir() / "test_init.mcool";
  const std::array<std::uint32_t, 4> resolutions{
      // clang-format off
      base_resolution,
      base_resolution * 2,
      base_resolution * 4,
      base_resolution * 8
      // clang-format on
  };

  SECTION("coarsen on construction") {
    SECTION("valid resolutions") {
      std::ignore = MultiResFile::create(path.string(), File(base_path.string()),
                                         resolutions.begin(), resolutions.end(), true);

      CHECK(utils::is_multires_file(path.string()));
    }

    SECTION("invalid resolutions") {
      std::vector<std::uint32_t> resolutions_{base_resolution / 2};
      CHECK_THROWS(MultiResFile::create(path.string(), File(base_path.string()),
                                        resolutions_.begin(), resolutions_.end(), true));
      resolutions_ = {base_resolution + 1};
      CHECK_THROWS(MultiResFile::create(path.string(), File(base_path.string()),
                                        resolutions_.begin(), resolutions_.end(), true));
    }
  }

  SECTION("construct then initialize") {
    const Reference chroms{Chromosome{0, "chr1", 10000}, Chromosome{1, "chr2", 5000}};
    auto mclr = MultiResFile::create(path.string(), chroms, true);

    for (const auto res : resolutions) {
      std::ignore = File::create(fmt::format(FMT_STRING("{}::/resolutions/{}"), path.string(), res),
                                 chroms, res);
    }

    CHECK(utils::is_multires_file(path.string()));
  }
}

TEST_CASE("MultiResCooler: create resolutions", "[cooler][short]") {
  const auto base_path = datadir / "cooler" / "cooler_test_file.cool";
  const File base_clr(base_path.string());
  const auto base_resolution = base_clr.resolution();

  const auto path = testdir() / "test_create_resolutions.mcool";
  const std::array<std::uint32_t, 3> resolutions{
      // clang-format off
      base_resolution * 2,
      base_resolution * 4,
      base_resolution * 8
      // clang-format on
  };

  auto mclr = MultiResFile::create(path.string(), base_clr.chromosomes(), true);
  mclr.copy_resolution(base_clr);

  SECTION("valid resolutions") {
    for (const auto res : resolutions) {
      std::ignore = mclr.create_resolution(res);
    }

    CHECK(mclr.resolutions().size() == resolutions.size() + 1);
  }
  SECTION("invalid resolutions") {
    CHECK_THROWS(mclr.create_resolution(base_resolution / 2));
    CHECK_THROWS(mclr.create_resolution(base_resolution + 1));
  }
}

TEST_CASE("MultiResCooler: normalizations", "[cooler][short]") {
  const auto path = datadir / "cooler" / "multires_cooler_test_file.mcool";

  const MultiResFile mclr{path.string()};

  namespace balancing = hictk::balancing;

  const std::vector expected_norms_union{balancing::Method{"weight1"},
                                         balancing::Method{"weight2"}};
  const std::vector expected_norms_intersection{balancing::Method{"weight1"}};

  auto found_norms = mclr.avail_normalizations("union");
  REQUIRE(expected_norms_union.size() == found_norms.size());
  for (std::size_t i = 0; i < expected_norms_union.size(); ++i) {
    CHECK(expected_norms_union[i] == found_norms[i]);
  }

  found_norms = mclr.avail_normalizations("intersection");
  REQUIRE(expected_norms_intersection.size() == found_norms.size());
  for (std::size_t i = 0; i < expected_norms_intersection.size(); ++i) {
    CHECK(expected_norms_intersection[i] == found_norms[i]);
  }

  CHECK_THROWS_AS(mclr.avail_normalizations("invalid"), std::invalid_argument);

  CHECK(&mclr.avail_normalizations("intersection") == &mclr.avail_normalizations("intersection"));
  CHECK(&mclr.avail_normalizations("intersection") == &mclr.avail_normalizations("union"));
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::cooler::test::multires_cooler_file
