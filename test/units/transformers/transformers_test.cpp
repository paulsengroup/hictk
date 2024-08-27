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

#ifdef HICTK_WITH_ARROW

namespace internal {

template <typename N>
[[nodiscard]] N get_scalar(const std::shared_ptr<arrow::ChunkedArray>& col, std::int64_t i) {
  assert(!!col);
  auto res = col->GetScalar(i);
  if (!res.ok()) {
    throw std::runtime_error(res.status().message());
  }

  if constexpr (std::is_same_v<N, std::string>) {
    res = std::static_pointer_cast<arrow::DictionaryScalar>(*res)->GetEncodedValue();
    if (!res.ok()) {
      throw std::runtime_error(res.status().message());
    }
    return (*res)->ToString();
  }

  if constexpr (std::is_same_v<N, std::uint32_t>) {
    return std::static_pointer_cast<arrow::UInt32Scalar>(*res)->value;
  }
  if constexpr (std::is_same_v<N, std::uint64_t>) {
    return std::static_pointer_cast<arrow::UInt64Scalar>(*res)->value;
  }
  if constexpr (std::is_same_v<N, std::int32_t>) {
    return std::static_pointer_cast<arrow::Int32Scalar>(*res)->value;
  }
  if constexpr (std::is_same_v<N, std::int64_t>) {
    return std::static_pointer_cast<arrow::Int64Scalar>(*res)->value;
  }

  if constexpr (std::is_floating_point_v<N>) {
    return std::static_pointer_cast<arrow::DoubleScalar>(*res)->value;
  }

  throw std::logic_error("not implemented");
}
}  // namespace internal

template <std::int64_t i, typename N>
static void compare_pixel(const std::shared_ptr<arrow::Table>& table, const ThinPixel<N>& p) {
  assert(!!table);

  REQUIRE(i < table->num_rows());

  CHECK(internal::get_scalar<std::uint64_t>(table->GetColumnByName("bin1_id"), i) == p.bin1_id);
  CHECK(internal::get_scalar<std::uint64_t>(table->GetColumnByName("bin2_id"), i) == p.bin2_id);
  CHECK(internal::get_scalar<N>(table->GetColumnByName("count"), i) == p.count);
}

template <std::int64_t i, typename N>  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void compare_pixel(const std::shared_ptr<arrow::Table>& table, const Pixel<N>& p) {
  assert(!!table);

  REQUIRE(i < table->num_rows());

  CHECK(internal::get_scalar<std::string>(table->GetColumnByName("chrom1"), i) ==
        p.coords.bin1.chrom().name());
  CHECK(internal::get_scalar<std::uint32_t>(table->GetColumnByName("start1"), i) ==
        p.coords.bin1.start());
  CHECK(internal::get_scalar<std::uint32_t>(table->GetColumnByName("end1"), i) ==
        p.coords.bin1.end());
  CHECK(internal::get_scalar<std::string>(table->GetColumnByName("chrom2"), i) ==
        p.coords.bin2.chrom().name());
  CHECK(internal::get_scalar<std::uint32_t>(table->GetColumnByName("start2"), i) ==
        p.coords.bin2.start());
  CHECK(internal::get_scalar<std::uint32_t>(table->GetColumnByName("end2"), i) ==
        p.coords.bin2.end());
  CHECK(internal::get_scalar<N>(table->GetColumnByName("count"), i) == p.count);
}

namespace internal {
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
[[nodiscard]] static std::vector<ThinPixel<std::uint8_t>> arrow_table_to_coo_vector(
    const std::shared_ptr<arrow::Table>& data) {
  assert(!!data);

  std::vector<ThinPixel<std::uint8_t>> buff(static_cast<std::size_t>(data->num_rows()));

  const auto bin1_ids = data->GetColumnByName("bin1_id");
  const auto bin2_ids = data->GetColumnByName("bin2_id");

  for (std::int64_t i = 0; i < data->num_rows(); ++i) {
    buff[static_cast<std::size_t>(i)].bin1_id = get_scalar<std::uint64_t>(bin1_ids, i);
    buff[static_cast<std::size_t>(i)].bin2_id = get_scalar<std::uint64_t>(bin2_ids, i);
  }
  return buff;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
[[nodiscard]] static std::vector<Pixel<std::uint8_t>> arrow_table_to_bg2_vector(
    const Reference& chroms, const std::shared_ptr<arrow::Table>& data) {
  assert(!!data);

  std::vector<Pixel<std::uint8_t>> buff(static_cast<std::size_t>(data->num_rows()));

  const auto chrom1_ids = data->GetColumnByName("chrom1");
  const auto start1 = data->GetColumnByName("start1");
  const auto end1 = data->GetColumnByName("end1");

  const auto chrom2_ids = data->GetColumnByName("chrom2");
  const auto start2 = data->GetColumnByName("start2");
  const auto end2 = data->GetColumnByName("end2");

  for (std::int64_t i = 0; i < data->num_rows(); ++i) {
    buff[static_cast<std::size_t>(i)] = Pixel{chroms.at(get_scalar<std::string>(chrom1_ids, i)),
                                              get_scalar<std::uint32_t>(start1, i),
                                              get_scalar<std::uint32_t>(end1, i),
                                              chroms.at(get_scalar<std::string>(chrom2_ids, i)),
                                              get_scalar<std::uint32_t>(start2, i),
                                              get_scalar<std::uint32_t>(end2, i),
                                              std::uint8_t{}};
  }
  return buff;
}
}  // namespace internal

template <DataFrameFormat format, QuerySpan span>
// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void validate_format(const Reference& chroms, const std::shared_ptr<arrow::Table>& table) {
  if constexpr (format == DataFrameFormat::COO) {
    const auto pixels = internal::arrow_table_to_coo_vector(table);
    if constexpr (span == QuerySpan::upper_triangle) {
      for (const auto& pixel : pixels) {
        CHECK(pixel.bin1_id <= pixel.bin2_id);
      }
    } else if constexpr (span == QuerySpan::lower_triangle) {
      for (const auto& pixel : pixels) {
        CHECK(pixel.bin1_id >= pixel.bin2_id);
      }
    } else {
      throw std::logic_error("not implemented");
    }
    return;
  }

  assert(format == DataFrameFormat::BG2);
  const auto pixels = internal::arrow_table_to_bg2_vector(chroms, table);
  if constexpr (span == QuerySpan::upper_triangle) {
    for (const auto& pixel : pixels) {
      CHECK(pixel.coords.bin1 <= pixel.coords.bin2);
    }
  } else if constexpr (span == QuerySpan::lower_triangle) {
    for (const auto& pixel : pixels) {
      CHECK(pixel.coords.bin1 >= pixel.coords.bin2);
    }
  } else {
    throw std::logic_error("not implemented");
  }
}

#endif

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
      using N = std::int32_t;
      const auto path = datadir / "cooler/ENCFF993FGR.2500000.cool";
      const cooler::File clr(path.string());
      const auto& bins = clr.bins();
      auto sel = clr.fetch("chr1");
      auto first = sel.begin<N>();
      auto last = sel.end<N>();

      SECTION("COO<int> upper_triangle") {
        constexpr auto format = DataFrameFormat::COO;
        constexpr auto span = QuerySpan::upper_triangle;
        const auto table = ToDataFrame(first, last, format, nullptr, span)();

        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(2)->type() == *arrow::int32());

        // check head
        compare_pixel<0>(table, ThinPixel<N>{0, 0, 266106});
        compare_pixel<1>(table, ThinPixel<N>{0, 1, 32868});
        compare_pixel<2>(table, ThinPixel<N>{0, 2, 13241});

        // check tail
        compare_pixel<4462>(table, ThinPixel<N>{98, 98, 1001844});
        compare_pixel<4463>(table, ThinPixel<N>{98, 99, 68621});
        compare_pixel<4464>(table, ThinPixel<N>{99, 99, 571144});

        validate_format<format, span>(clr.chromosomes(), table);
      }

      SECTION("COO<int> lower_triangle") {
        constexpr auto format = DataFrameFormat::COO;
        constexpr auto span = QuerySpan::lower_triangle;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, clr.bins_ptr(), span)();

        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(2)->type() == *arrow::int32());

        // check head
        compare_pixel<0>(table, ThinPixel<N>{0, 0, 266106});
        compare_pixel<1>(table, ThinPixel<N>{1, 0, 32868});
        compare_pixel<2>(table, ThinPixel<N>{1, 1, 375662});

        // check tail
        compare_pixel<4462>(table, ThinPixel<N>{99, 97, 24112});
        compare_pixel<4463>(table, ThinPixel<N>{99, 98, 68621});
        compare_pixel<4464>(table, ThinPixel<N>{99, 99, 571144});

        validate_format<format, span>(clr.chromosomes(), table);
      }

      SECTION("COO<int> full") {
        constexpr auto format = DataFrameFormat::COO;
        constexpr auto span = QuerySpan::full;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, clr.bins_ptr(), span)();

        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 8'836);
        CHECK(*table->column(2)->type() == *arrow::int32());

        // check head
        compare_pixel<0>(table, ThinPixel<N>{0, 0, 266106});
        compare_pixel<1>(table, ThinPixel<N>{0, 1, 32868});
        compare_pixel<2>(table, ThinPixel<N>{0, 2, 13241});

        // check tail
        compare_pixel<8833>(table, ThinPixel<N>{99, 97, 24112});
        compare_pixel<8834>(table, ThinPixel<N>{99, 98, 68621});
        compare_pixel<8835>(table, ThinPixel<N>{99, 99, 571144});
      }

      SECTION("BG2<int> upper_triangle") {
        constexpr auto format = DataFrameFormat::BG2;
        constexpr auto span = QuerySpan::upper_triangle;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, clr.bins_ptr(), span)();

        CHECK(table->num_columns() == 7);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(6)->type() == *arrow::int32());

        // check head
        compare_pixel<0>(table, Pixel<N>{bins.at("chr1", 0), bins.at("chr1", 0), 266106});
        compare_pixel<1>(table, Pixel<N>{bins.at("chr1", 0), bins.at("chr1", 2'500'000), 32868});
        compare_pixel<2>(table, Pixel<N>{bins.at("chr1", 0), bins.at("chr1", 5'000'000), 13241});

        // check tail
        compare_pixel<4462>(
            table, Pixel<N>{bins.at("chr1", 245'000'000), bins.at("chr1", 245'000'000), 1001844});
        compare_pixel<4463>(
            table, Pixel<N>{bins.at("chr1", 245'000'000), bins.at("chr1", 247'500'000), 68621});
        compare_pixel<4464>(
            table, Pixel<N>{bins.at("chr1", 247'500'000), bins.at("chr1", 247'500'000), 571144});

        validate_format<format, span>(clr.chromosomes(), table);
      }

      SECTION("BG2<int> lower_triangle") {
        constexpr auto format = DataFrameFormat::BG2;
        constexpr auto span = QuerySpan::lower_triangle;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, clr.bins_ptr(), span)();

        CHECK(table->num_columns() == 7);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(6)->type() == *arrow::int32());

        // check head
        compare_pixel<0>(table, Pixel<N>{bins.at("chr1", 0), bins.at("chr1", 0), 266106});
        compare_pixel<1>(table, Pixel<N>{bins.at("chr1", 2'500'000), bins.at("chr1", 0), 32868});
        compare_pixel<2>(table,
                         Pixel<N>{bins.at("chr1", 2'500'000), bins.at("chr1", 2'500'000), 375662});

        // check tail
        compare_pixel<4462>(
            table, Pixel<N>{bins.at("chr1", 247'500'000), bins.at("chr1", 242'500'000), 24112});
        compare_pixel<4463>(
            table, Pixel<N>{bins.at("chr1", 247'500'000), bins.at("chr1", 245'000'000), 68621});
        compare_pixel<4464>(
            table, Pixel<N>{bins.at("chr1", 247'500'000), bins.at("chr1", 247'500'000), 571144});

        validate_format<format, span>(clr.chromosomes(), table);
      }

      SECTION("BG2<int> full") {
        constexpr auto format = DataFrameFormat::BG2;
        constexpr auto span = QuerySpan::full;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, clr.bins_ptr(), span)();

        CHECK(table->num_columns() == 7);
        CHECK(table->num_rows() == 8'836);
        CHECK(*table->column(6)->type() == *arrow::int32());

        // check head
        compare_pixel<0>(table, Pixel<N>{bins.at("chr1", 0), bins.at("chr1", 0), 266106});
        compare_pixel<1>(table, Pixel<N>{bins.at("chr1", 0), bins.at("chr1", 2'500'000), 32868});
        compare_pixel<2>(table, Pixel<N>{bins.at("chr1", 0), bins.at("chr1", 5'000'000), 13241});

        // check tail
        compare_pixel<8833>(
            table, Pixel<N>{bins.at("chr1", 247'500'000), bins.at("chr1", 242'500'000), 24112});
        compare_pixel<8834>(
            table, Pixel<N>{bins.at("chr1", 247'500'000), bins.at("chr1", 245'000'000), 68621});
        compare_pixel<8835>(
            table, Pixel<N>{bins.at("chr1", 247'500'000), bins.at("chr1", 247'500'000), 571144});
      }

      SECTION("COO<float> upper_triangle") {
        constexpr auto format = DataFrameFormat::COO;
        constexpr auto span = QuerySpan::upper_triangle;

        auto first_fp = sel.begin<double>();
        auto last_fp = sel.end<double>();
        const auto table = ToDataFrame(first_fp, last_fp, format, nullptr, span)();

        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 4'465);
        CHECK(*table->column(2)->type() == *arrow::float64());

        // check head
        compare_pixel<0>(table, ThinPixel<double>{0, 0, 266106.0});
        compare_pixel<1>(table, ThinPixel<double>{0, 1, 32868.0});
        compare_pixel<2>(table, ThinPixel<double>{0, 2, 13241.0});

        // check tail
        compare_pixel<4462>(table, ThinPixel<double>{98, 98, 1001844.0});
        compare_pixel<4463>(table, ThinPixel<double>{98, 99, 68621.0});
        compare_pixel<4464>(table, ThinPixel<double>{99, 99, 571144.0});

        validate_format<format, span>(clr.chromosomes(), table);
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

  if constexpr (TEST_TO_DATAFRAME) {
    SECTION("ToDataFrame") {
      using N = std::int32_t;
      const hic::File hf(path.string(), 2'500'000);
      auto sel = hf.fetch("chr2L");
      auto first = sel.begin<N>();
      auto last = sel.end<N>();

      SECTION("COO<int> upper_triangle") {
        constexpr auto format = DataFrameFormat::COO;
        constexpr auto span = QuerySpan::upper_triangle;
        const auto table = ToDataFrame(first, last, format, nullptr, span)();

        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 55);
        CHECK(*table->column(2)->type() == *arrow::int32());

        validate_format<format, span>(hf.chromosomes(), table);
      }

      SECTION("COO<int> lower_triangle") {
        constexpr auto format = DataFrameFormat::COO;
        constexpr auto span = QuerySpan::lower_triangle;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, hf.bins_ptr(), span)();

        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 55);
        CHECK(*table->column(2)->type() == *arrow::int32());

        validate_format<format, span>(hf.chromosomes(), table);
      }

      SECTION("COO<int> full") {
        constexpr auto format = DataFrameFormat::COO;
        constexpr auto span = QuerySpan::full;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, hf.bins_ptr(), span)();

        CHECK(table->num_columns() == 3);
        CHECK(table->num_rows() == 100);
        CHECK(*table->column(2)->type() == *arrow::int32());
      }

      SECTION("BG2<int> upper_triangle") {
        constexpr auto format = DataFrameFormat::BG2;
        constexpr auto span = QuerySpan::upper_triangle;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, hf.bins_ptr(), span)();

        CHECK(table->num_columns() == 7);
        CHECK(table->num_rows() == 55);
        CHECK(*table->column(6)->type() == *arrow::int32());

        validate_format<format, span>(hf.chromosomes(), table);
      }

      SECTION("BG2<int> lower_triangle") {
        constexpr auto format = DataFrameFormat::BG2;
        constexpr auto span = QuerySpan::lower_triangle;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, hf.bins_ptr(), span)();

        CHECK(table->num_columns() == 7);
        CHECK(table->num_rows() == 55);
        CHECK(*table->column(6)->type() == *arrow::int32());

        validate_format<format, span>(hf.chromosomes(), table);
      }

      SECTION("BG2<int> full") {
        constexpr auto format = DataFrameFormat::BG2;
        constexpr auto span = QuerySpan::full;
        CHECK_THROWS(ToDataFrame(first, last, format, nullptr, span));
        const auto table = ToDataFrame(first, last, format, hf.bins_ptr(), span)();

        CHECK(table->num_columns() == 7);
        CHECK(table->num_rows() == 100);
        CHECK(*table->column(6)->type() == *arrow::int32());
      }
      SECTION("empty range") {
        const auto table = ToDataFrame(last, last)();
        CHECK(table->num_rows() == 0);
      }
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
