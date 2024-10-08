// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <filesystem>
#include <vector>

#include "hictk/balancing/weights.hpp"
#include "hictk/cooler/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

TEST_CASE("Cooler: read weights", "[cooler][short]") {
  const auto path1 = datadir / "cooler_test_file.cool";
  const auto path2 = datadir / "ENCFF993FGR.2500000.cool";

  const cooler::File clr1{path1.string()};
  const cooler::File clr2{path2.string()};

  SECTION("wo/ weights") { CHECK(clr1.avail_normalizations().empty()); }
  SECTION("w/ weights") {
    CHECK(clr2.avail_normalizations().size() == 8);
    CHECK(clr2.has_normalization("SCALE"));
    CHECK(!clr2.has_normalization("FOOBAR"));

    CHECK(clr2.normalization("SCALE").type() == hictk::balancing::Weights::Type::DIVISIVE);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: write weights", "[cooler][short]") {
  const auto path1 = datadir / "cooler_test_file.cool";
  auto path2 = testdir() / "cooler_test_write_weights1.cool";
  auto path3 = testdir() / "cooler_test_write_weights2.cool";

  std::filesystem::remove(path2);
  std::filesystem::remove(path3);
  std::filesystem::copy(path1, path2);
  REQUIRE_FALSE(File(path2.string()).has_normalization("weight"));

  const auto num_bins = File(path1.string()).bins().size();

  SECTION("correct shape") {
    const std::vector<double> weights(num_bins, 1.23);
    File::write_weights(path2.string(), "weight", weights.begin(), weights.end());

    const auto w = File(path2.string()).normalization("weight");
    CHECK(w.size() == weights.size());
  }

  SECTION("incorrect shape") {
    std::vector<double> weights{};
    CHECK_THROWS(File::write_weights(path2.string(), "weight", weights.begin(), weights.end()));

    weights.resize(num_bins - 1);
    CHECK_THROWS(File::write_weights(path2.string(), "weight", weights.begin(), weights.end()));

    weights.resize(num_bins + 1);
    CHECK_THROWS(File::write_weights(path2.string(), "weight", weights.begin(), weights.end()));
  }

  SECTION("invalid name") {
    std::vector<double> weights{};
    CHECK_THROWS(File::write_weights(path2.string(), "", weights.begin(), weights.end()));
  }

  SECTION("overwriting") {
    const std::vector<double> weights(num_bins, 1.23);
    File::write_weights(path2.string(), "weight", weights.begin(), weights.end());

    File::write_weights(path2.string(), "weight", weights.begin(), weights.end(), true);

    CHECK_THROWS(
        File::write_weights(path2.string(), "weight", weights.begin(), weights.end(), false));
  }

  SECTION("write on file creation") {
    const File fin(path1.string());
    auto fout = File::create(path3.string(), fin.chromosomes(), fin.resolution());

    const std::vector<double> weights(num_bins, 1.23);
    fout.write_weights("weight", weights.begin(), weights.end());
    fout.write_weights("weight2", weights.begin(), weights.end());
  }

  SECTION("attempt write on read-only file") {
    constexpr std::array<double, 1> w{};
    CHECK_THROWS(File(path2.string()).write_weights("weights", w.begin(), w.end()));
  }
}

}  // namespace hictk::cooler::test::cooler_file
