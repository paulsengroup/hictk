// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/hic_file_stream.hpp"

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <filesystem>
#include <string>

using namespace hictk;

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data/hic"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

const auto pathV8 = test::datadir / "4DNFIZ1ZVXC8.hic8";
const auto pathV9 = test::datadir / "4DNFIZ1ZVXC8.hic9";

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("readHeader (v8)", "[hic][v8][short]") {
  constexpr std::array<std::int32_t, 10> resolutions{2500000, 1000000, 500000, 250000, 100000,
                                                     50000,   25000,   10000,  5000,   1000};
  constexpr auto* genomeID = "dm6";
  constexpr auto nChromosomes = 9;

  const auto header = internal::HiCFileStream(pathV8).header();
  CHECK(header.url == pathV8);
  CHECK(header.masterIndexOffset == 131515430);
  CHECK(header.genomeID == genomeID);
  CHECK(header.nChromosomes() == nChromosomes);
  CHECK(header.version == 8);
  CHECK(header.nviPosition == -1);
  CHECK(header.nviLength == -1);

  REQUIRE(header.nResolutions() == resolutions.size());
  for (std::size_t i = 0; i < resolutions.size(); ++i) {
    CHECK(resolutions[i] == header.resolutions[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("readHeader (v9)", "[hic][v9][short]") {
  constexpr std::array<std::int32_t, 10> resolutions{2500000, 1000000, 500000, 250000, 100000,
                                                     50000,   25000,   10000,  5000,   1000};
  constexpr auto* genomeID = "dm6";
  constexpr auto nChromosomes = 9;

  const auto header = internal::HiCFileStream(pathV9).header();

  CHECK(header.url == pathV9);
  CHECK(header.masterIndexOffset == 130706734);
  CHECK(header.genomeID == genomeID);
  CHECK(header.nChromosomes() == nChromosomes);
  CHECK(header.version == 9);
  CHECK(header.nviPosition == 131417220);
  CHECK(header.nviLength == 6600);

  REQUIRE(header.nResolutions() == resolutions.size());
  for (std::size_t i = 0; i < resolutions.size(); ++i) {
    CHECK(resolutions[i] == header.resolutions[i]);
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("readFooter (v8)", "[hic][v8][short]") {
  internal::HiCFileStream s(pathV8);
  const auto chr2L = s.header().chromosomes.at("chr2L");
  const auto chr2R = s.header().chromosomes.at("chr2R");
  // first 5 expected values
  constexpr std::array<double, 5> expected1{864.6735714977542, 620.9907283534235, 311.1254999778368,
                                            203.9822974509631, 147.9273228359822};
  // last 5 expected values
  constexpr std::array<double, 5> expected2{0.008417076032024847, 0.008417076032024847,
                                            0.008417076032024847, 0.008417076032024847,
                                            0.008417076032024847};

  SECTION("observed NONE BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2L.index, MatrixType::observed,
                                NormalizationMethod::NONE, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::NONE);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340697);
    CHECK(f.c1Norm().empty());
    CHECK(f.c2Norm().empty());
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed VC BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2R.index, MatrixType::observed,
                                NormalizationMethod::VC, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::VC);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11389664);
    CHECK(f.c1Norm().size() == 4703);
    CHECK(f.c2Norm().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed VC_SQRT BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2R.index, MatrixType::observed,
                                NormalizationMethod::VC_SQRT, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::VC_SQRT);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11389664);
    CHECK(f.c1Norm().size() == 4703);
    CHECK(f.c2Norm().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed KR BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2R.index, MatrixType::observed,
                                NormalizationMethod::KR, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::KR);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11389664);
    CHECK(f.c1Norm().size() == 4703);
    CHECK(f.c2Norm().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed SCALE BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2R.index, MatrixType::observed,
                                NormalizationMethod::SCALE, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::SCALE);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11389664);
    CHECK(f.c1Norm().size() == 4703);
    CHECK(f.c2Norm().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("oe NONE BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2L.index, MatrixType::oe, NormalizationMethod::NONE,
                                MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::oe);
    CHECK(f.normalization() == NormalizationMethod::NONE);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340697);
    CHECK(f.c1Norm().empty());
    CHECK(f.c2Norm().empty());
    REQUIRE(f.expectedValues().size() == 6415);

    for (std::size_t i = 0; i < expected1.size(); ++i) {
      const auto j = f.expectedValues().size() - (expected2.size() - i);
      CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues()[i]));
      CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues()[j]));
    }
  }

  SECTION("expected NONE BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2L.index, MatrixType::expected,
                                NormalizationMethod::NONE, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::expected);
    CHECK(f.normalization() == NormalizationMethod::NONE);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340697);
    CHECK(f.c1Norm().empty());
    CHECK(f.c2Norm().empty());
    REQUIRE(f.expectedValues().size() == 6415);

    for (std::size_t i = 0; i < expected1.size(); ++i) {
      const auto j = f.expectedValues().size() - (expected2.size() - i);
      CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues()[i]));
      CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues()[j]));
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("readFooter (v9)", "[hic][v9][short]") {
  internal::HiCFileStream s(pathV9);
  const auto chr2L = s.header().chromosomes.at("chr2L");
  const auto chr2R = s.header().chromosomes.at("chr2R");
  // first 5 expected values
  constexpr std::array<double, 5> expected1{864.6735708339686, 620.990715491172, 311.1255023627755,
                                            203.9822882714327, 147.9273192507429};
  // last 5 expected values
  constexpr std::array<double, 5> expected2{0.008417075820557469, 0.008417075820557469,
                                            0.008417075820557469, 0.008417075820557469,
                                            0.008417075820557469};

  SECTION("observed NONE BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2L.index, MatrixType::observed,
                                NormalizationMethod::NONE, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::NONE);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340696);
    CHECK(f.c1Norm().empty());
    CHECK(f.c2Norm().empty());
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed VC BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2R.index, MatrixType::observed,
                                NormalizationMethod::VC, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::VC);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11625116);
    CHECK(f.c1Norm().size() == 4703);
    CHECK(f.c2Norm().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("observed VC_SQRT BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2R.index, MatrixType::observed,
                                NormalizationMethod::VC_SQRT, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::VC_SQRT);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11625116);
    CHECK(f.c1Norm().size() == 4703);
    CHECK(f.c2Norm().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  /*  TODO: for some reason KR normalization is missing
  SECTION("observed KR BP 5000") {
      const auto f = s.readFooter(chr2L.index, chr2R.index, MatrixType::observed,
  NormalizationMethod::KR, MatrixUnit::BP, 5000);

      CHECK(f.matrix_type() == MatrixType::observed);
      CHECK(f.normalization() == NormalizationMethod::KR);
      CHECK(f.unit() == MatrixUnit::BP);
      CHECK(f.resolution() == 5000);
      CHECK(f.fileOffset() == 11625116);
      CHECK(f.c1Norm().size() == 4703);
      CHECK(f.c2Norm().size() == 5058);
      CHECK(f._expectedValues().empty());
  } */

  SECTION("observed SCALE BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2R.index, MatrixType::observed,
                                NormalizationMethod::SCALE, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::observed);
    CHECK(f.normalization() == NormalizationMethod::SCALE);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 11625116);
    CHECK(f.c1Norm().size() == 4703);
    CHECK(f.c2Norm().size() == 5058);
    CHECK(f.expectedValues().empty());
  }

  SECTION("oe NONE BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2L.index, MatrixType::oe, NormalizationMethod::NONE,
                                MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::oe);
    CHECK(f.normalization() == NormalizationMethod::NONE);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340696);
    CHECK(f.c1Norm().empty());
    CHECK(f.c2Norm().empty());
    REQUIRE(f.expectedValues().size() == 6415);

    for (std::size_t i = 0; i < expected1.size(); ++i) {
      const auto j = f.expectedValues().size() - (expected2.size() - i);
      CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues()[i]));
      CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues()[j]));
    }
  }

  SECTION("expected NONE BP 5000") {
    const auto f = s.readFooter(chr2L.index, chr2L.index, MatrixType::expected,
                                NormalizationMethod::NONE, MatrixUnit::BP, 5000);

    CHECK(f.matrix_type() == MatrixType::expected);
    CHECK(f.normalization() == NormalizationMethod::NONE);
    CHECK(f.unit() == MatrixUnit::BP);
    CHECK(f.resolution() == 5000);
    CHECK(f.fileOffset() == 340696);
    CHECK(f.c1Norm().empty());
    CHECK(f.c2Norm().empty());
    REQUIRE(f.expectedValues().size() == 6415);

    for (std::size_t i = 0; i < expected1.size(); ++i) {
      const auto j = f.expectedValues().size() - (expected2.size() - i);
      CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues()[i]));
      CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues()[j]));
    }
  }
}
