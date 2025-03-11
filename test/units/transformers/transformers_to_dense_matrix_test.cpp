// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef HICTK_WITH_EIGEN

#include <algorithm>
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
#include "hictk/transformers/to_dense_matrix.hpp"

namespace hictk::test::transformers {

using namespace hictk::transformers;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
template <typename Matrix>
[[nodiscard]] static double sum_finite(const Matrix& matrix) noexcept {
  const double* first = matrix.data();
  const double* last = first + matrix.size();  // NOLINT(*-pointer-arithmetic)
  return std::accumulate(first, last, 0.0, [](double accumulator, double n) {
    return accumulator + (std::isfinite(n) ? n : 0.0);
  });
}

template <typename Matrix>
[[nodiscard]] static std::size_t count_nans(const Matrix& matrix) noexcept {
  // NOLINTNEXTLINE(*-pointer-arithmetic)
  return static_cast<std::size_t>(std::count_if(matrix.data(), matrix.data() + matrix.size(),
                                                [&](auto n) { return std::isnan(n); }));
}

TEST_CASE("Transformers (cooler): to dense matrix", "[transformers][short]") {
  SECTION("cis full") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch("chr1"), std::int32_t{}, QuerySpan::full)();
    CHECK(matrix.rows() == 100);
    CHECK(matrix.cols() == 100);
    CHECK(matrix.sum() == 140'900'545);
    CHECK(matrix == matrix.transpose());
  }

  SECTION("cis upper_triangle") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix =
        ToDenseMatrix(clr.fetch("chr1"), std::int32_t{}, QuerySpan::upper_triangle)();
    CHECK(matrix.rows() == 100);
    CHECK(matrix.cols() == 100);
    CHECK(matrix.sum() == 112'660'799);
    CHECK(matrix.isUpperTriangular());
  }

  SECTION("cis lower_triangle") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix =
        ToDenseMatrix(clr.fetch("chr1"), std::int32_t{}, QuerySpan::lower_triangle)();
    CHECK(matrix.rows() == 100);
    CHECK(matrix.cols() == 100);
    CHECK(matrix.sum() == 112'660'799);
    CHECK(matrix.isLowerTriangular());
  }

  SECTION("(cis, asymmetric) full") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix =
        ToDenseMatrix(clr.fetch("chr1:192,565,354-202,647,735", "chr1:197,313,124-210,385,543"),
                      std::int32_t{}, QuerySpan::full)();
    CHECK(matrix.rows() == 5);
    CHECK(matrix.cols() == 7);
    CHECK(matrix.sum() == 5'426'501);
  }

  SECTION("(cis, normalized) full") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix =
        ToDenseMatrix(clr.fetch("chr1", balancing::Method{"VC"}), 0.0, QuerySpan::full)();
    CHECK(matrix.rows() == 100);
    CHECK(matrix.cols() == 100);

    CHECK_THAT(sum_finite(matrix), Catch::Matchers::WithinRel(140900543.1839076));
    CHECK(count_nans(matrix) == 1164);
  }

  SECTION("trans upper_triangle") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix =
        ToDenseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{}, QuerySpan::upper_triangle)();
    CHECK(matrix.rows() == 100);
    CHECK(matrix.cols() == 97);
    CHECK(matrix.sum() == 6'413'076);
  }

  SECTION("trans lower_triangle") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    CHECK_THROWS(
        ToDenseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{}, QuerySpan::lower_triangle));
  }

  SECTION("trans full") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{}, QuerySpan::full)();
    CHECK(matrix.rows() == 100);
    CHECK(matrix.cols() == 97);
    CHECK(matrix.sum() == 6'413'076);
  }

  SECTION("(trans, normalized) full") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix =
        ToDenseMatrix(clr.fetch("chr1", "chr2", balancing::Method{"VC"}), 0.0, QuerySpan::full)();
    CHECK(matrix.rows() == 100);
    CHECK(matrix.cols() == 97);

    CHECK_THAT(sum_finite(matrix), Catch::Matchers::WithinRel(6185975.980057132));
    CHECK(count_nans(matrix) == 582);
  }

  SECTION("gw full") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch(), std::uint32_t{}, QuerySpan::full)();
    CHECK(matrix.rows() == 1249);
    CHECK(matrix.cols() == 1249);
    CHECK(matrix.sum() == 2'671'244'699);
  }

  SECTION("gw upper_triangle") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch(), std::int32_t{}, QuerySpan::upper_triangle)();
    CHECK(matrix.rows() == 1249);
    CHECK(matrix.cols() == 1249);
    CHECK(matrix.sum() == 1'868'866'491);
    CHECK(matrix.isUpperTriangular());
  }

  SECTION("gw lower_triangle") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch(), std::int32_t{}, QuerySpan::lower_triangle)();
    CHECK(matrix.rows() == 1249);
    CHECK(matrix.cols() == 1249);
    CHECK(matrix.sum() == 1'868'866'491);
    CHECK(matrix.isLowerTriangular());
  }

  SECTION("gw full (storage-mode=square)") {
    const auto path =
        datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/1000";
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch(), std::uint32_t{}, QuerySpan::full)();
    CHECK(matrix.rows() == 3000);
    CHECK(matrix.cols() == 3000);
    CHECK(matrix.sum() == 594'006'205);
  }

  SECTION("gw upper_triangle (storage-mode=square)") {
    const auto path =
        datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/1000";
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch(), std::int32_t{}, QuerySpan::upper_triangle)();
    CHECK(matrix.rows() == 3000);
    CHECK(matrix.cols() == 3000);
    CHECK(matrix.sum() == 336'795'259);
    CHECK(matrix.isUpperTriangular());
  }

  SECTION("gw lower_triangle (storage-mode=square)") {
    const auto path =
        datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/1000";
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch(), std::int32_t{}, QuerySpan::lower_triangle)();
    CHECK(matrix.rows() == 3000);
    CHECK(matrix.cols() == 3000);
    CHECK(matrix.sum() == 257'471'326);
    CHECK(matrix.isLowerTriangular());
  }

  SECTION("gw w/ diagonal band") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    constexpr std::uint64_t diagonal_band_width{10};
    const cooler::File clr(path.string());
    const auto matrix = ToDenseMatrix(clr.fetch(balancing::Method::NONE(), true), std::uint32_t{},
                                      QuerySpan::full, diagonal_band_width)();
    CHECK(matrix.rows() == 1249);
    CHECK(matrix.cols() == 1249);
    CHECK(matrix.sum() == 1'539'111'295);
  }

  SECTION("regression PR #154") {
    const auto path = datadir / "cooler" / "cooler_test_file.cool";
    const cooler::File clr(path.string());
    const auto matrix =
        ToDenseMatrix(clr.fetch("1:0-5,000,000", "1:2,500,000-7,500,000"), std::int32_t{})();

    CHECK(matrix.rows() == 50);
    CHECK(matrix.cols() == 50);
    CHECK(matrix.sum() == 442);
  }

  SECTION("invalid queries") {
    const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());

    CHECK_THROWS(ToDenseMatrix(clr.fetch("chr1", "chr2"), 0, QuerySpan::lower_triangle));
    CHECK_THROWS(ToDenseMatrix(clr.fetch("chr1", balancing::Method{"weight"}), 0));
  }
}

TEST_CASE("Transformers (hic): to dense matrix", "[transformers][short]") {
  const auto path = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();

  SECTION("cis") {
    const hic::File hf(path, 2'500'000);
    const auto matrix = ToDenseMatrix(hf.fetch("chr2L"), std::int32_t{})();
    CHECK(matrix.rows() == 10);
    CHECK(matrix.cols() == 10);
    CHECK(matrix.sum() == 22'929'541);
    CHECK(matrix == matrix.transpose());
  }

  SECTION("cis, normalized") {
    const hic::File hf(path, 2'500'000);
    const auto matrix = ToDenseMatrix(hf.fetch("chr2L", balancing::Method{"VC"}), 0.0)();
    CHECK(matrix.rows() == 10);
    CHECK(matrix.cols() == 10);

    CHECK_THAT(sum_finite(matrix), Catch::Matchers::WithinRel(22929540.99999999, 1.0e-6));
    CHECK(count_nans(matrix) == 0);
  }

  SECTION("trans") {
    const hic::File hf(path, 2'500'000);
    const auto matrix = ToDenseMatrix(hf.fetch("chr2L", "chr2R"), std::int32_t{})();
    CHECK(matrix.rows() == 10);
    CHECK(matrix.cols() == 11);
    CHECK(matrix.sum() == 1'483'112);
  }

  SECTION("gw") {
    const hic::File hf(path, 2'500'000);
    const auto matrix = ToDenseMatrix(hf.fetch(), std::int32_t{})();
    CHECK(matrix.rows() == 60);
    CHECK(matrix.cols() == 60);
    CHECK(matrix.sum() == 149'078'427);
    CHECK(matrix == matrix.transpose());
  }

  SECTION("gw, normalized") {
    const hic::File hf(path, 2'500'000);
    const auto matrix = ToDenseMatrix(hf.fetch(balancing::Method{"VC"}), 0.0)();
    CHECK(matrix.rows() == 60);
    CHECK(matrix.cols() == 60);

    CHECK_THAT(sum_finite(matrix), Catch::Matchers::WithinRel(146874129.31714758, 1.0e-6));
    CHECK(count_nans(matrix) == 119);
  }

  SECTION("invalid queries") {
    const hic::File hf(path, 2'500'000);

    CHECK_THROWS(ToDenseMatrix(hf.fetch("chr2L", "chr2R"), 0, QuerySpan::lower_triangle));
    CHECK_THROWS(ToDenseMatrix(hf.fetch("chr2L", balancing::Method::VC()), 0));
  }
}

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::transformers

#endif
