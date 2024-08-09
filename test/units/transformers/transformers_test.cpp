// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <parallel_hashmap/btree.h>

#include <array>
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
#include "hictk/transformers/coarsen.hpp"
#include "hictk/transformers/join_genomic_coords.hpp"
#include "hictk/transformers/pixel_merger.hpp"
#include "hictk/transformers/stats.hpp"
#include "hictk/transformers/to_dataframe.hpp"
#include "hictk/transformers/to_dense_matrix.hpp"
#include "hictk/transformers/to_sparse_matrix.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

#ifdef HICTK_WITH_ARROW
constexpr auto TEST_TO_DATAFRAME = true;
#else
constexpr auto TEST_TO_DATAFRAME = false;
#endif

#ifdef HICTK_WITH_EIGEN
constexpr auto TEST_TO_SPARSE_MATRIX = true;
constexpr auto TEST_TO_DENSE_MATRIX = true;
#else
constexpr auto TEST_TO_SPARSE_MATRIX = false;
constexpr auto TEST_TO_DENSE_MATRIX = false;
#endif

using namespace hictk::transformers;

struct Coords {
  std::uint64_t bin1{};  // NOLINT
  std::uint64_t bin2{};  // NOLINT

  bool operator==(const Coords& other) const noexcept {
    return bin1 == other.bin1 && bin2 == other.bin2;
  }
  bool operator<(const Coords& other) const noexcept {
    if (bin1 == other.bin1) {
      return bin2 < other.bin2;
    }
    return bin1 < other.bin1;
  }
};

template <typename PixelIt>
static phmap::btree_map<Coords, std::int32_t> merge_pixels_hashmap(
    const std::vector<PixelIt>& heads, const std::vector<PixelIt>& tails) {
  assert(heads.size() == tails.size());

  phmap::btree_map<Coords, std::int32_t> map{};
  for (std::size_t i = 0; i < heads.size(); ++i) {
    std::for_each(heads[i], tails[i], [&](const ThinPixel<std::int32_t>& p) {
      const Coords c{p.bin1_id, p.bin2_id};

      if (map.contains(c)) {
        map[c] += p.count;
      } else {
        map[c] = p.count;
      }
    });
  }
  return map;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (cooler)", "[transformers][short]") {
  SECTION("join genomic coords") {
    const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());

    auto sel = clr.fetch("chr1", 5'000'000, 10'000'000);
    SECTION("range with data") {
      const auto jsel =
          JoinGenomicCoords(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr.bins_ptr());
      constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
      const auto pixels = jsel.read_all();
      REQUIRE(pixels.size() == expected.size());
      for (std::size_t i = 0; i < expected.size(); ++i) {
        CHECK(pixels[i].coords.bin1.start() == expected[i]);
      }
    }
    SECTION("empty range") {
      const auto jsel =
          JoinGenomicCoords(sel.end<std::int32_t>(), sel.end<std::int32_t>(), clr.bins_ptr());
      CHECK(jsel.begin() == jsel.end());
    }
  }

  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
  SECTION("PixelMerger", "[transformers][short]") {
    const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";

    const cooler::File clr(path.string());
    const auto sel1 = clr.fetch("chr1:0-100,000,000");
    const auto sel2 = clr.fetch("chr1:50,000,000-150,000,000");
    const auto sel3 = clr.fetch("chr2:50,000,000-150,000,000");

    using It = decltype(sel1.template begin<std::int32_t>());

    SECTION("range with data") {
      const std::vector<It> heads{sel1.template begin<std::int32_t>(),
                                  sel2.template begin<std::int32_t>(),
                                  sel3.template begin<std::int32_t>()};
      const std::vector<It> tails{sel1.template end<std::int32_t>(),
                                  sel2.template end<std::int32_t>(),
                                  sel3.template end<std::int32_t>()};

      const transformers::PixelMerger<It> merger(heads, tails);
      const auto pixels = transformers::PixelMerger<It>(heads, tails).read_all();
      const auto expected_pixels = merge_pixels_hashmap(heads, tails);

      REQUIRE(pixels.size() == expected_pixels.size());

      for (const auto& p : pixels) {
        REQUIRE(expected_pixels.contains({p.bin1_id, p.bin2_id}));
        CHECK(expected_pixels.at({p.bin1_id, p.bin2_id}) == p.count);
      }
    }

    SECTION("one iterator") {
      const std::vector<It> heads{sel1.template begin<std::int32_t>()};
      const std::vector<It> tails{sel1.template end<std::int32_t>()};
      const transformers::PixelMerger<It> merger(heads, tails);
      const auto pixels = transformers::PixelMerger<It>(heads, tails).read_all();
      const auto expected_pixels = merge_pixels_hashmap(heads, tails);

      REQUIRE(pixels.size() == expected_pixels.size());
      for (const auto& p : pixels) {
        REQUIRE(expected_pixels.contains({p.bin1_id, p.bin2_id}));
        CHECK(expected_pixels.at({p.bin1_id, p.bin2_id}) == p.count);
      }
    }

    SECTION("empty range") {
      const std::vector<It> heads{sel1.template begin<std::int32_t>(),
                                  sel2.template end<std::int32_t>(),
                                  sel3.template begin<std::int32_t>()};
      const std::vector<It> tails{sel1.template end<std::int32_t>(),
                                  sel2.template end<std::int32_t>(),
                                  sel3.template end<std::int32_t>()};
      const transformers::PixelMerger<It> merger(heads, tails);
      const auto pixels = transformers::PixelMerger<It>(heads, tails).read_all();
      const auto expected_pixels = merge_pixels_hashmap(heads, tails);

      REQUIRE(pixels.size() == expected_pixels.size());
      for (const auto& p : pixels) {
        REQUIRE(expected_pixels.contains({p.bin1_id, p.bin2_id}));
        CHECK(expected_pixels.at({p.bin1_id, p.bin2_id}) == p.count);
      }
    }

    SECTION("no iterators") {
      const std::vector<It> heads{};
      const std::vector<It> tails{};
      const transformers::PixelMerger<It> merger(heads, tails);

      CHECK(merger.begin() == merger.end());
    }
  }

  SECTION("coarsen") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/200000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch("1");
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = clr2.fetch("1");

    const auto v1 = sel1.read_all();
    const auto v2 = sel2.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }
  SECTION("coarsen recursive") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/400000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch("1");
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = CoarsenPixels(sel1.begin(), sel1.end(), sel1.dest_bins_ptr(), 2);
    auto sel3 = clr2.fetch("1");

    const auto v1 = sel2.read_all();
    const auto v2 = sel3.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  SECTION("coarsen gw") {
    const auto path1 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const auto path2 = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/200000";
    const cooler::File clr1(path1.string());
    const cooler::File clr2(path2.string());

    auto sel = clr1.fetch();
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);
    auto sel2 = clr2.fetch();

    const auto v1 = sel1.read_all();
    const auto v2 = sel2.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  SECTION("coarsen empty range") {
    const auto path = datadir / "cooler/multires_cooler_test_file.mcool::/resolutions/100000";
    const cooler::File clr1(path.string());

    auto sel = clr1.fetch();
    auto sel1 = CoarsenPixels(sel.end<std::int32_t>(), sel.end<std::int32_t>(), clr1.bins_ptr(), 2);

    CHECK(sel1.begin() == sel1.end());
  }

  SECTION("stats") {
    const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    auto sel = clr.fetch("chr1");
    auto first = sel.begin<std::int32_t>();
    auto last = sel.end<std::int32_t>();

    SECTION("range with data") {
      CHECK_THAT(avg(first, last), Catch::Matchers::WithinRel(25231.981858902574));
      CHECK(nnz(first, last) == 4'465);
      CHECK(max(first, last) == 1'357'124);
      CHECK(sum(first, last) == 112'660'799);
    }

    SECTION("empty range") {
      CHECK(avg(last, last) == 0);
      CHECK(nnz(last, last) == 0);
      CHECK(max(last, last) == 0);
      CHECK(sum(last, last) == 0);
    }
  }

  if constexpr (TEST_TO_DATAFRAME) {
    SECTION("ToDataFrame") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      auto sel = clr.fetch("chr1");
      auto first = sel.begin<std::int32_t>();
      auto last = sel.end<std::int32_t>();

      auto get_int_scalar = [](const auto& column, const auto i) {
        return std::static_pointer_cast<arrow::Int32Scalar>(column->GetScalar(i).MoveValueUnsafe())
            ->value;
      };

      auto get_float_scalar = [](const auto& column, const auto i) {
        return std::static_pointer_cast<arrow::DoubleScalar>(column->GetScalar(i).MoveValueUnsafe())
            ->value;
      };

      SECTION("COO<int> wo/ transpose") {
        const auto table = ToDataFrame(first, last, DataFrameFormat::COO, nullptr, false)();
        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(2)->type() == *arrow::int32());

        // check head
        CHECK(get_int_scalar(table->column(2), 0) == 266106);
        CHECK(get_int_scalar(table->column(2), 1) == 32868);
        CHECK(get_int_scalar(table->column(2), 2) == 13241);

        // check tail
        CHECK(get_int_scalar(table->column(2), 4462) == 1001844);
        CHECK(get_int_scalar(table->column(2), 4463) == 68621);
        CHECK(get_int_scalar(table->column(2), 4464) == 571144);
      }

      SECTION("COO<int> w/ transpose") {
        const auto table = ToDataFrame(first, last, DataFrameFormat::COO, nullptr, true)();
        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(2)->type() == *arrow::int32());

        // check head
        CHECK(get_int_scalar(table->column(2), 0) == 266106);
        CHECK(get_int_scalar(table->column(2), 1) == 32868);
        CHECK(get_int_scalar(table->column(2), 2) == 375662);

        // check tail
        CHECK(get_int_scalar(table->column(2), 4462) == 24112);
        CHECK(get_int_scalar(table->column(2), 4463) == 68621);
        CHECK(get_int_scalar(table->column(2), 4464) == 571144);
      }

      SECTION("BG2<int> wo/ transpose") {
        const auto table = ToDataFrame(first, last, DataFrameFormat::BG2, clr.bins_ptr(), false)();
        CHECK(table->num_columns() == 7);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(6)->type() == *arrow::int32());

        CHECK(get_int_scalar(table->column(6), 0) == 266106);
        CHECK(get_int_scalar(table->column(6), 1) == 32868);
        CHECK(get_int_scalar(table->column(6), 2) == 13241);

        // check tail
        CHECK(get_int_scalar(table->column(6), 4462) == 1001844);
        CHECK(get_int_scalar(table->column(6), 4463) == 68621);
        CHECK(get_int_scalar(table->column(6), 4464) == 571144);
      }

      SECTION("BG2<int> w/ transpose") {
        const auto table = ToDataFrame(first, last, DataFrameFormat::BG2, clr.bins_ptr(), true)();
        CHECK(table->num_columns() == 7);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(6)->type() == *arrow::int32());

        // check head
        CHECK(get_int_scalar(table->column(6), 0) == 266106);
        CHECK(get_int_scalar(table->column(6), 1) == 32868);
        CHECK(get_int_scalar(table->column(6), 2) == 375662);

        // check tail
        CHECK(get_int_scalar(table->column(6), 4462) == 24112);
        CHECK(get_int_scalar(table->column(6), 4463) == 68621);
        CHECK(get_int_scalar(table->column(6), 4464) == 571144);
      }

      SECTION("COO<float> wo/ transpose") {
        auto first_fp = sel.begin<double>();
        auto last_fp = sel.end<double>();
        const auto table = ToDataFrame(first_fp, last_fp, DataFrameFormat::COO, nullptr, false)();
        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(2)->type() == *arrow::float64());

        // check head
        CHECK(get_float_scalar(table->column(2), 0) == 266106.0);
        CHECK(get_float_scalar(table->column(2), 1) == 32868.0);
        CHECK(get_float_scalar(table->column(2), 2) == 13241.0);

        // check tail
        CHECK(get_float_scalar(table->column(2), 4462) == 1001844.0);
        CHECK(get_float_scalar(table->column(2), 4463) == 68621.0);
        CHECK(get_float_scalar(table->column(2), 4464) == 571144.0);
      }

      SECTION("empty range") {
        const auto table = ToDataFrame(last, last)();
        CHECK(table->num_rows() == 0);
      }
    }
  }

  if constexpr (TEST_TO_SPARSE_MATRIX) {
    SECTION("ToSparseMatrix (cis) wo/ transpose") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1"), std::int32_t{}, false)();
      CHECK(matrix.nonZeros() == 4465);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 112'660'799);
      CHECK(matrix.triangularView<Eigen::StrictlyLower>().sum() == 0);
    }

    SECTION("ToSparseMatrix (cis) w/ transpose") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1"), std::int32_t{}, true)();
      CHECK(matrix.nonZeros() == 4465);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 112'660'799);
      CHECK(matrix.triangularView<Eigen::StrictlyUpper>().sum() == 0);
    }

    SECTION("ToSparseMatrix (trans) wo/ transpose") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{}, false)();
      CHECK(matrix.nonZeros() == 9118);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 97);
      CHECK(matrix.sum() == 6'413'076);
    }

    SECTION("ToSparseMatrix (trans) w/ transpose") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{}, true)();
      CHECK(matrix.nonZeros() == 9118);
      CHECK(matrix.rows() == 97);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 6'413'076);
    }

    SECTION("ToSparseMatrix (gw) wo/ transpose") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch(), std::int32_t{}, false)();
      CHECK(matrix.nonZeros() == 718'781);
      CHECK(matrix.rows() == 1249);
      CHECK(matrix.cols() == 1249);
      CHECK(matrix.sum() == 1'868'866'491);
      CHECK(matrix.triangularView<Eigen::StrictlyLower>().sum() == 0);
    }

    SECTION("ToSparseMatrix (gw) w/ transpose") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch(), std::int32_t{}, true)();
      CHECK(matrix.nonZeros() == 718'781);
      CHECK(matrix.rows() == 1249);
      CHECK(matrix.cols() == 1249);
      CHECK(matrix.sum() == 1'868'866'491);
      CHECK(matrix.triangularView<Eigen::StrictlyUpper>().sum() == 0);
    }
  }

  if constexpr (TEST_TO_DENSE_MATRIX) {
    SECTION("ToDenseMatrix (cis) w/ mirroring") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToDenseMatrix(clr.fetch("chr1"), std::int32_t{}, true)();
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 140'900'545);
      CHECK(matrix == matrix.transpose());
    }

    SECTION("ToDenseMatrix (cis) wo/ mirroring") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToDenseMatrix(clr.fetch("chr1"), std::int32_t{}, false)();
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 112'660'799);
    }

    SECTION("ToDenseMatrix (trans)") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToDenseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{})();
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 97);
      CHECK(matrix.sum() == 6'413'076);
    }

    SECTION("ToDenseMatrix (gw) w/ mirroring") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToDenseMatrix(clr.fetch(), std::uint32_t{}, true)();
      CHECK(matrix.rows() == 1249);
      CHECK(matrix.cols() == 1249);
      CHECK(matrix.sum() == 2'671'244'699);
    }

    SECTION("ToDenseMatrix (gw) wo/ mirroring") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToDenseMatrix(clr.fetch(), std::int32_t{}, false)();
      CHECK(matrix.rows() == 1249);
      CHECK(matrix.cols() == 1249);
      CHECK(matrix.sum() == 1'868'866'491);
    }

    SECTION("ToDenseMatrix regression PR #154") {
      const auto path = datadir / "cooler/cooler_test_file.cool";
      const cooler::File clr(path.string());
      const auto matrix =
          ToDenseMatrix(clr.fetch("1:0-5,000,000", "1:2,500,000-7,500,000"), std::int32_t{})();

      CHECK(matrix.rows() == 50);
      CHECK(matrix.cols() == 50);
      CHECK(matrix.sum() == 442);
    }
  }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("Transformers (hic)", "[transformers][short]") {
  auto path = datadir / "hic/4DNFIZ1ZVXC8.hic8";

  SECTION("join genomic coords") {
    const hic::File hf(path.string(), 2'500'000);
    auto sel = hf.fetch("chr2L", 5'000'000, 10'000'000);
    const auto jsel =
        JoinGenomicCoords(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), hf.bins_ptr());
    constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
    const auto pixels = jsel.read_all();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(pixels[i].coords.bin1.start() == expected[i]);
    }
  }

  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
  SECTION("PixelMerger", "[transformers][short]") {
    const hic::File hf(path.string(), 100'000);
    const auto sel1 = hf.fetch("chr2L:0-10,000,000");
    const auto sel2 = hf.fetch("chr2L:5,000,000-15,000,000");
    const auto sel3 = hf.fetch("chr2R:5,000,000-15,000,000");

    using It = decltype(sel1.template begin<std::int32_t>());
    const std::vector<It> heads{sel1.template begin<std::int32_t>(),
                                sel2.template begin<std::int32_t>(),
                                sel3.template begin<std::int32_t>()};
    const std::vector<It> tails{sel1.template end<std::int32_t>(),
                                sel2.template end<std::int32_t>(),
                                sel3.template end<std::int32_t>()};

    const transformers::PixelMerger<It> merger(heads, tails);
    const auto pixels = transformers::PixelMerger<It>(heads, tails).read_all();
    const auto expected_pixels = merge_pixels_hashmap(heads, tails);

    REQUIRE(pixels.size() == expected_pixels.size());

    for (const auto& p : pixels) {
      REQUIRE(expected_pixels.contains({p.bin1_id, p.bin2_id}));
      CHECK(expected_pixels.at({p.bin1_id, p.bin2_id}) == p.count);
    }
  }

  SECTION("coarsen") {
    const hic::File hf1(path.string(), 500'000);
    const hic::File hf2(path.string(), 2'500'000);

    auto sel = hf1.fetch("chr2R");
    auto sel1 =
        CoarsenPixels(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), hf1.bins_ptr(), 5);
    auto sel2 = hf2.fetch("chr2R");

    const auto v1 = sel1.read_all();
    const auto v2 = sel2.read_all<std::int32_t>();
    REQUIRE(v1.size() == v2.size());

    for (std::size_t i = 0; i < v1.size(); ++i) {
      CHECK(v1[i] == v2[i].to_thin());
    }
  }

  if constexpr (TEST_TO_SPARSE_MATRIX) {
    SECTION("ToSparseMatrix (cis)") {
      const hic::File hf(path.string(), 2'500'000);
      const auto matrix = ToSparseMatrix(hf.fetch("chr2L"), std::int32_t{})();
      CHECK(matrix.nonZeros() == 55);
      CHECK(matrix.rows() == 10);
      CHECK(matrix.cols() == 10);
      CHECK(matrix.sum() == 19'968'156);
    }

    SECTION("ToSparseMatrix (trans)") {
      const hic::File hf(path.string(), 2'500'000);
      const auto matrix = ToSparseMatrix(hf.fetch("chr2L", "chr2R"), std::int32_t{})();
      CHECK(matrix.nonZeros() == 110);
      CHECK(matrix.rows() == 10);
      CHECK(matrix.cols() == 11);
      CHECK(matrix.sum() == 1'483'112);
    }

    SECTION("ToSparseMatrix (gw)") {
      const hic::File hf(path.string(), 2'500'000);
      const auto matrix = ToSparseMatrix(hf.fetch(), std::int32_t{})();
      CHECK(matrix.nonZeros() == 1770);
      CHECK(matrix.rows() == 60);
      CHECK(matrix.cols() == 60);
      CHECK(matrix.sum() == 119'208'613);
    }
  }

  if constexpr (TEST_TO_DENSE_MATRIX) {
    SECTION("ToDenseMatrix (cis)") {
      const hic::File hf(path.string(), 2'500'000);
      const auto matrix = ToDenseMatrix(hf.fetch("chr2L"), std::int32_t{})();
      CHECK(matrix.rows() == 10);
      CHECK(matrix.cols() == 10);
      CHECK(matrix.sum() == 22'929'541);
      CHECK(matrix == matrix.transpose());
    }

    SECTION("ToDenseMatrix (trans)") {
      const hic::File hf(path.string(), 2'500'000);
      const auto matrix = ToDenseMatrix(hf.fetch("chr2L", "chr2R"), std::int32_t{})();
      CHECK(matrix.rows() == 10);
      CHECK(matrix.cols() == 11);
      CHECK(matrix.sum() == 1'483'112);
    }

    SECTION("ToDenseMatrix (gw)") {
      const hic::File hf(path.string(), 2'500'000);
      const auto matrix = ToDenseMatrix(hf.fetch(), std::int32_t{})();
      CHECK(matrix.rows() == 60);
      CHECK(matrix.cols() == 60);
      CHECK(matrix.sum() == 149'078'427);
      CHECK(matrix == matrix.transpose());
    }
  }
}

}  // namespace hictk::test::transformers
