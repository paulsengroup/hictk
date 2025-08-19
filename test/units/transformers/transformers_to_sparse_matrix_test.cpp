// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef HICTK_WITH_EIGEN

#include <spdlog/spdlog.h>

#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/pixel.hpp"
#include "hictk/test/testdir.hpp"
#include "hictk/transformers/to_sparse_matrix.hpp"

namespace hictk::test::transformers {

using namespace hictk::transformers;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
TEST_CASE("Transformers (cooler): to sparse matrix", "[transformers][short]") {
  spdlog::set_level(spdlog::level::trace);

  for (const auto low_mem : {true, false}) {
    const auto* suffix = low_mem ? "low-mem" : "fast";

    SECTION(fmt::format(FMT_STRING("cis upper_triangle ({})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch("chr1"), std::int32_t{}, QuerySpan::upper_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 4465);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 112'660'799);
      CHECK(matrix.triangularView<Eigen::StrictlyLower>().sum() == 0);
    }

    SECTION(fmt::format(FMT_STRING("cis lower_triangle ({})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch("chr1"), std::int32_t{}, QuerySpan::lower_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 4465);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 112'660'799);
      CHECK(matrix.triangularView<Eigen::StrictlyUpper>().sum() == 0);
    }

    SECTION(fmt::format(FMT_STRING("cis full ({})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch("chr1"), std::int32_t{}, QuerySpan::full, low_mem)();
      CHECK(matrix.nonZeros() == 8836);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 140'900'545);
      CHECK(matrix.triangularView<Eigen::Upper>().sum() ==
            matrix.triangularView<Eigen::Lower>().sum());
    }

    SECTION(fmt::format(FMT_STRING("cis upper_triangle (asymmetric; {})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1:0-10,000,000", "chr1:0-21,000,000"),
                                         std::int32_t{}, QuerySpan::upper_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 30);
      CHECK(matrix.rows() == 4);
      CHECK(matrix.cols() == 9);
      CHECK(matrix.sum() == 2'231'517);
    }

    SECTION(fmt::format(FMT_STRING("cis lower_triangle (asymmetric; {})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1:0-10,000,000", "chr1:0-21,000,000"),
                                         std::int32_t{}, QuerySpan::lower_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 10);
      CHECK(matrix.rows() == 4);
      CHECK(matrix.cols() == 9);
      CHECK(matrix.sum() == 2'007'400);
    }

    SECTION(fmt::format(FMT_STRING("cis full (asymmetric; {})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1:0-10,000,000", "chr1:0-21,000,000"),
                                         std::int32_t{}, QuerySpan::full, low_mem)();
      CHECK(matrix.nonZeros() == 36);
      CHECK(matrix.rows() == 4);
      CHECK(matrix.cols() == 9);
      CHECK(matrix.sum() == 2'411'797);
    }

    SECTION(fmt::format(FMT_STRING("trans upper_triangle ({})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{},
                                         QuerySpan::upper_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 9118);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 97);
      CHECK(matrix.sum() == 6'413'076);
    }

    SECTION(fmt::format(FMT_STRING("trans lower_triangle ({})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      CHECK_THROWS(ToSparseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{},
                                  QuerySpan::lower_triangle, low_mem));
    }

    SECTION(fmt::format(FMT_STRING("trans full ({})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{}, QuerySpan::full, low_mem)();
      CHECK(matrix.nonZeros() == 9118);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 97);
      CHECK(matrix.sum() == 6'413'076);
    }

    SECTION(fmt::format(FMT_STRING("gw upper_triangle ({})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch(), std::int32_t{}, QuerySpan::upper_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 718'781);
      CHECK(matrix.rows() == 1249);
      CHECK(matrix.cols() == 1249);
      CHECK(matrix.sum() == 1'868'866'491);
      CHECK(matrix.triangularView<Eigen::StrictlyLower>().sum() == 0);
    }

    SECTION(fmt::format(FMT_STRING("gw lower_triangle ({})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch(), std::int32_t{}, QuerySpan::lower_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 718'781);
      CHECK(matrix.rows() == 1249);
      CHECK(matrix.cols() == 1249);
      CHECK(matrix.sum() == 1'868'866'491);
      CHECK(matrix.triangularView<Eigen::StrictlyUpper>().sum() == 0);
    }

    SECTION(fmt::format(FMT_STRING("gw full (storage-mode=square; {})"), suffix)) {
      const auto path =
          datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/1000";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch(), std::uint32_t{}, QuerySpan::full, low_mem)();
      CHECK(matrix.nonZeros() == 4'241'909);
      CHECK(matrix.rows() == 3000);
      CHECK(matrix.cols() == 3000);
      CHECK(matrix.sum() == 594'006'205);
    }

    SECTION(fmt::format(FMT_STRING("gw upper_triangle (storage-mode=square; {})"), suffix)) {
      const auto path =
          datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/1000";
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch(), std::int32_t{}, QuerySpan::upper_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 2'423'572);
      CHECK(matrix.rows() == 3000);
      CHECK(matrix.cols() == 3000);
      CHECK(matrix.sum() == 336'795'259);
      CHECK(matrix.triangularView<Eigen::StrictlyLower>().sum() == 0);
    }

    SECTION(fmt::format(FMT_STRING("gw lower_triangle (storage-mode=square; {})"), suffix)) {
      const auto path =
          datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/1000";
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch(), std::int32_t{}, QuerySpan::lower_triangle, low_mem)();
      CHECK(matrix.nonZeros() == 1'820'117);
      CHECK(matrix.rows() == 3000);
      CHECK(matrix.cols() == 3000);
      CHECK(matrix.sum() == 257'471'326);
      CHECK(matrix.triangularView<Eigen::StrictlyUpper>().sum() == 0);
    }

    SECTION(fmt::format(FMT_STRING("gw full (diagonal band=10; {})"), suffix)) {
      const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
      constexpr std::uint64_t diagonal_band_width{10};
      const cooler::File clr(path.string());
      const auto matrix =
          ToSparseMatrix(clr.fetch(balancing::Method::NONE(), true), std::uint32_t{},
                         QuerySpan::full, low_mem, diagonal_band_width)();
      CHECK(matrix.rows() == 1249);
      CHECK(matrix.cols() == 1249);
      CHECK(matrix.sum() == 1'539'111'295);
    }
  }
  SECTION("invalid queries") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());

    CHECK_THROWS(ToSparseMatrix(clr.fetch("chr1", "chr2"), 0, QuerySpan::lower_triangle));
    CHECK_THROWS(ToSparseMatrix(clr.fetch("chr1", balancing::Method::VC()), 0));
  }
}

TEST_CASE("Transformers (hic): to sprase matrix", "[transformers][short]") {
  const auto path = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();

  SECTION("cis") {
    const hic::File hf(path, 2'500'000);
    const auto matrix = ToSparseMatrix(hf.fetch("chr2L"), std::int32_t{})();
    CHECK(matrix.nonZeros() == 55);
    CHECK(matrix.rows() == 10);
    CHECK(matrix.cols() == 10);
    CHECK(matrix.sum() == 19'968'156);
  }

  SECTION("trans") {
    const hic::File hf(path, 2'500'000);
    const auto matrix = ToSparseMatrix(hf.fetch("chr2L", "chr2R"), std::int32_t{})();
    CHECK(matrix.nonZeros() == 110);
    CHECK(matrix.rows() == 10);
    CHECK(matrix.cols() == 11);
    CHECK(matrix.sum() == 1'483'112);
  }

  SECTION("gw") {
    const hic::File hf(path, 2'500'000);
    const auto matrix = ToSparseMatrix(hf.fetch(), std::int32_t{})();
    CHECK(matrix.nonZeros() == 1770);
    CHECK(matrix.rows() == 60);
    CHECK(matrix.cols() == 60);
    CHECK(matrix.sum() == 119'208'613);
  }

  SECTION("invalid queries") {
    const hic::File hf(path, 2'500'000);

    CHECK_THROWS(ToSparseMatrix(hf.fetch("chr2L", "chr2R"), 0, QuerySpan::lower_triangle));
    CHECK_THROWS(ToSparseMatrix(hf.fetch("chr2L", balancing::Method::VC()), 0));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::transformers

#endif
