// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <filesystem>

#include "hictk/cooler.hpp"
#include "tmpdir.hpp"

namespace hictk::cooler::test::cooler_file {

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Cooler: write weights", "[cooler][short]") {
  auto path1 = datadir / "cooler_test_file.cool";
  auto path2 = testdir() / "cooler_test_write_weights1.cool";
  auto path3 = testdir() / "cooler_test_write_weights2.cool";

  std::filesystem::remove(path2);
  std::filesystem::remove(path3);
  std::filesystem::copy(path1, path2);
  REQUIRE_FALSE(File::open_read_only(path2.string()).has_weights("weight"));

  const auto num_bins = File::open_read_only(path1.string()).bins().size();

  SECTION("correct shape") {
    const std::vector<double> weights(num_bins, 1.23);
    File::write_weights(path2.string(), "weight", weights.begin(), weights.end());

    const auto w = *File::open_read_only(path2.string()).read_weights("weight");
    CHECK(w().size() == weights.size());
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
    const auto fin = File::open_read_only(path1.string());
    auto fout = File::create_new_cooler(path3.string(), fin.chromosomes(), fin.bin_size());

    const std::vector<double> weights(num_bins, 1.23);
    fout.write_weights("weight", weights.begin(), weights.end());
    fout.write_weights("weight2", weights.begin(), weights.end());
  }

  SECTION("attempt write on read-only file") {
    constexpr std::array<double, 1> w{};
    CHECK_THROWS(File::open_read_only(path2.string()).write_weights("weights", w.begin(), w.end()));
  }
}

}  // namespace hictk::cooler::test::cooler_file
