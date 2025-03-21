// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/file_reader.hpp"

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/test/testdir.hpp"

using namespace hictk::hic;

namespace hictk::hic::test::file_reader {

static const auto& datadir = hictk::test::datadir;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)

// NOLINTNEXTLINE(cert-err58-cpp)
const auto pathV8 = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();
// NOLINTNEXTLINE(cert-err58-cpp)
const auto pathV9 = (datadir / "hic" / "4DNFIZ1ZVXC8.hic9").string();

[[maybe_unused]] static void check_weights_are_constant(const balancing::Weights& weights,
                                                        double value = 1.0) {
  REQUIRE(weights.is_constant());
  REQUIRE(!weights.empty());
  CHECK_THAT(*weights.begin(), Catch::Matchers::WithinRel(value));
}

TEST_CASE("HiC: read header (v8)", "[hic][v8][short]") {
  constexpr std::array<std::uint32_t, 10> resolutions{
      1'000, 5'000, 10'000, 25'000, 50'000, 100'000, 250'000, 500'000, 1'000'000, 2'500'000};
  constexpr auto* genomeID = "dm6";
  constexpr auto nChromosomes = 9;

  const auto header = internal::HiCFileReader(pathV8).header();
  CHECK(header.url == pathV8);
  CHECK(header.footerPosition == 131515430);
  CHECK(header.genomeID == genomeID);
  CHECK(header.chromosomes.size() == nChromosomes);
  CHECK(header.version == 8);
  CHECK(header.normVectorIndexPosition == -1);
  CHECK(header.normVectorIndexLength == -1);

  REQUIRE(header.resolutions.size() == resolutions.size());
  for (std::size_t i = 0; i < resolutions.size(); ++i) {
    CHECK(resolutions[i] == header.resolutions[i]);
  }
}

TEST_CASE("HiC: read header (v9)", "[hic][v9][short]") {
  constexpr std::array<std::uint32_t, 10> resolutions{
      1'000, 5'000, 10'000, 25'000, 50'000, 100'000, 250'000, 500'000, 1'000'000, 2'500'000};
  constexpr auto* genomeID = "dm6";
  constexpr auto nChromosomes = 9;

  const auto header = internal::HiCFileReader(pathV9).header();

  CHECK(header.url == pathV9);
  CHECK(header.footerPosition == 130706734);
  CHECK(header.genomeID == genomeID);
  CHECK(header.chromosomes.size() == nChromosomes);
  CHECK(header.version == 9);
  CHECK(header.normVectorIndexPosition == 131417220);
  CHECK(header.normVectorIndexLength == 6600);

  REQUIRE(header.resolutions.size() == resolutions.size());
  for (std::size_t i = 0; i < resolutions.size(); ++i) {
    CHECK(resolutions[i] == header.resolutions[i]);
  }
}

TEST_CASE("HiC: read footer (v8)", "[hic][v8][short]") {
  internal::HiCFileReader s(pathV8);
  const auto chr2L = s.header().chromosomes.at("chr2L");
  const auto chr2R = s.header().chromosomes.at("chr2R");
  const BinTable bins(s.header().chromosomes, 5000);
  // first 5 expected values
  constexpr std::array<double, 5> expected1{864.6735714977542, 620.9907283534235, 311.1254999778368,
                                            203.9822974509631, 147.9273228359822};
  // last 5 expected values
  constexpr std::array<double, 5> expected2{0.008417076032024847, 0.008417076032024847,
                                            0.008417076032024847, 0.008417076032024847,
                                            0.008417076032024847};

  SECTION("observed NONE BP 5000") {
    auto w = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2L, bins, MatrixType::observed,
                                 hictk::balancing::Method::NONE(), MatrixUnit::BP, w, w);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == "NONE");
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340697);
    check_weights_are_constant(f.weights1());
    check_weights_are_constant(f.weights2());
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed VC BP 5000") {
    auto w1 = std::make_shared<balancing::Weights>();
    auto w2 = std::make_shared<balancing::Weights>();

    const auto f = s.read_footer(chr2L, chr2R, bins, MatrixType::observed,
                                 hictk::balancing::Method::VC(), MatrixUnit::BP, w1, w2);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == "VC");
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11389664);
    CHECK(f.weights1().size() == 4703);
    CHECK(f.weights2().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed VC_SQRT BP 5000") {
    auto w1 = std::make_shared<balancing::Weights>();
    auto w2 = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2R, bins, MatrixType::observed,
                                 hictk::balancing::Method::VC_SQRT(), MatrixUnit::BP, w1, w2);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == hictk::balancing::Method::VC_SQRT());
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11389664);
    CHECK(f.weights1().size() == 4703);
    CHECK(f.weights2().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed KR BP 5000") {
    auto w1 = std::make_shared<balancing::Weights>();
    auto w2 = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2R, bins, MatrixType::observed,
                                 hictk::balancing::Method::KR(), MatrixUnit::BP, w1, w2);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == hictk::balancing::Method::KR());
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11389664);
    CHECK(f.weights1().size() == 4703);
    CHECK(f.weights2().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed SCALE BP 5000") {
    auto w1 = std::make_shared<balancing::Weights>();
    auto w2 = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2R, bins, MatrixType::observed,
                                 hictk::balancing::Method::SCALE(), MatrixUnit::BP, w1, w2);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == hictk::balancing::Method::SCALE());
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11389664);
    CHECK(f.weights1().size() == 4703);
    CHECK(f.weights2().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("oe NONE BP 5000") {
    auto w = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2L, bins, MatrixType::oe,
                                 hictk::balancing::Method::NONE(), MatrixUnit::BP, w, w);

    CHECK(f.matrix_type() == MatrixType::oe);
    CHECK(f.normalization() == "NONE");
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340697);
    check_weights_are_constant(f.weights1());
    check_weights_are_constant(f.weights2());
    REQUIRE(f.expectedValues().size() == 6415);

    for (std::size_t i = 0; i < expected1.size(); ++i) {
      const auto j = f.expectedValues().size() - (expected2.size() - i);
      CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues()[i]));
      CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues()[j]));
    }
  }

  SECTION("expected NONE BP 5000") {
    auto w = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2L, bins, MatrixType::expected,
                                 hictk::balancing::Method::NONE(), MatrixUnit::BP, w, w);

    CHECK(f.matrix_type() == MatrixType::expected);
    CHECK(f.normalization() == "NONE");
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340697);
    check_weights_are_constant(f.weights1());
    check_weights_are_constant(f.weights2());
    REQUIRE(f.expectedValues().size() == 6415);

    for (std::size_t i = 0; i < expected1.size(); ++i) {
      const auto j = f.expectedValues().size() - (expected2.size() - i);
      CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues()[i]));
      CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues()[j]));
    }
  }
}

TEST_CASE("HiC: read footer (v9)", "[hic][v9][short]") {
  internal::HiCFileReader s(pathV9);
  const auto chr2L = s.header().chromosomes.at("chr2L");
  const auto chr2R = s.header().chromosomes.at("chr2R");
  const BinTable bins(s.header().chromosomes, 5000);
  // first 5 expected values
  constexpr std::array<double, 5> expected1{864.6735708339686, 620.990715491172, 311.1255023627755,
                                            203.9822882714327, 147.9273192507429};
  // last 5 expected values
  constexpr std::array<double, 5> expected2{0.008417075820557469, 0.008417075820557469,
                                            0.008417075820557469, 0.008417075820557469,
                                            0.008417075820557469};

  SECTION("observed NONE BP 5000") {
    auto w = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2L, bins, MatrixType::observed,
                                 hictk::balancing::Method::NONE(), MatrixUnit::BP, w, w);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == "NONE");
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340696);
    check_weights_are_constant(f.weights1());
    check_weights_are_constant(f.weights2());
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed VC BP 5000") {
    auto w1 = std::make_shared<balancing::Weights>();
    auto w2 = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2R, bins, MatrixType::observed,
                                 hictk::balancing::Method::VC(), MatrixUnit::BP, w1, w2);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == hictk::balancing::Method::VC());
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11625116);
    CHECK(f.weights1().size() == 4703);
    CHECK(f.weights2().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed VC_SQRT BP 5000") {
    auto w1 = std::make_shared<balancing::Weights>();
    auto w2 = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2R, bins, MatrixType::observed,
                                 hictk::balancing::Method::VC_SQRT(), MatrixUnit::BP, w1, w2);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == hictk::balancing::Method::VC_SQRT());
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11625116);
    CHECK(f.weights1().size() == 4703);
    CHECK(f.weights2().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed SCALE BP 5000") {
    auto w1 = std::make_shared<balancing::Weights>();
    auto w2 = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2R, bins, MatrixType::observed,
                                 hictk::balancing::Method::SCALE(), MatrixUnit::BP, w1, w2);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == hictk::balancing::Method::SCALE());
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11625116);
    CHECK(f.weights1().size() == 4703);
    CHECK(f.weights2().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("oe NONE BP 5000") {
    auto w = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2L, bins, MatrixType::oe,
                                 hictk::balancing::Method::NONE(), MatrixUnit::BP, w, w);

    CHECK(f.matrix_type() == MatrixType::oe);
    CHECK(f.normalization() == "NONE");
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340696);
    check_weights_are_constant(f.weights1());
    check_weights_are_constant(f.weights2());
    REQUIRE(f.expectedValues().size() == 6415);

    for (std::size_t i = 0; i < expected1.size(); ++i) {
      const auto j = f.expectedValues().size() - (expected2.size() - i);
      CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues()[i]));
      CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues()[j]));
    }
  }

  SECTION("expected NONE BP 5000") {
    auto w = std::make_shared<balancing::Weights>();
    const auto f = s.read_footer(chr2L, chr2L, bins, MatrixType::expected,
                                 hictk::balancing::Method::NONE(), MatrixUnit::BP, w, w);

    CHECK(f.matrix_type() == MatrixType::expected);
    CHECK(f.normalization() == "NONE");
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340696);
    check_weights_are_constant(f.weights1());
    check_weights_are_constant(f.weights2());
    REQUIRE(f.expectedValues().size() == 6415);

    for (std::size_t i = 0; i < expected1.size(); ++i) {
      const auto j = f.expectedValues().size() - (expected2.size() - i);
      CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues()[i]));
      CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues()[j]));
    }
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::hic::test::file_reader
