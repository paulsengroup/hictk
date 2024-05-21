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
#include "hictk/transformers/to_sparse_matrix.hpp"
#include "hictk/transformers/to_dense_matrix.hpp"

namespace hictk::test {
inline const std::filesystem::path datadir{"test/data"};  // NOLINT(cert-err58-cpp)
}  // namespace hictk::test

namespace hictk::test::transformers {

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
    const auto jsel =
        JoinGenomicCoords(sel.begin<std::int32_t>(), sel.end<std::int32_t>(), clr.bins_ptr());
    constexpr std::array<std::uint32_t, 3> expected{5'000'000, 5'000'000, 7'500'000};
    const auto pixels = jsel.read_all();
    REQUIRE(pixels.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
      CHECK(pixels[i].coords.bin1.start() == expected[i]);
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

  SECTION("stats") {
    const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
    const cooler::File clr(path.string());
    auto sel = clr.fetch("chr1");
    auto first = sel.begin<std::int32_t>();
    auto last = sel.end<std::int32_t>();

    CHECK_THAT(avg(first, last), Catch::Matchers::WithinRel(25231.981858902574));
    CHECK(nnz(first, last) == 4'465);
    CHECK(max(first, last) == 1'357'124);
    CHECK(sum(first, last) == 112'660'799);
  }

  if constexpr (TEST_TO_SPARSE_MATRIX) {
    SECTION("ToSparseMatrix (cis)") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1"), std::int32_t{})();
      CHECK(matrix.nonZeros() == 4465);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 112'660'799);
    }

    SECTION("ToSparseMatrix (trans)") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToSparseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{})();
      CHECK(matrix.nonZeros() == 9118);
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 97);
      CHECK(matrix.sum() == 6'413'076);
    }
  }

  if constexpr (TEST_TO_DENSE_MATRIX) {
    SECTION("ToDenseMatrix (cis)") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToDenseMatrix(clr.fetch("chr1"), std::int32_t{})();
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 100);
      CHECK(matrix.sum() == 140'900'545);
    }

    SECTION("ToDenseMatrix (trans)") {
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto matrix = ToDenseMatrix(clr.fetch("chr1", "chr2"), std::int32_t{})();
      CHECK(matrix.rows() == 100);
      CHECK(matrix.cols() == 97);
      CHECK(matrix.sum() == 6'413'076);
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
    }
  }
}

}  // namespace hictk::test::transformers
