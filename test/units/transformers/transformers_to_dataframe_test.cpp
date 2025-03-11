// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef HICTK_WITH_ARROW

#include <arrow/scalar.h>
#include <arrow/table.h>

#include <cassert>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/pixel.hpp"
#include "hictk/test/testdir.hpp"
#include "hictk/transformers/to_dataframe.hpp"

namespace hictk::test::transformers {

using namespace hictk::transformers;

// NOLINTBEGIN(*-avoid-magic-numbers, readability-function-cognitive-complexity)
namespace internal {

template <typename N>
[[nodiscard]] static N get_scalar(const std::shared_ptr<arrow::ChunkedArray>& col, std::int64_t i) {
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

  CHECK(internal::get_scalar<std::int64_t>(table->GetColumnByName("bin1_id"), i) ==
        static_cast<std::int64_t>(p.bin1_id));
  CHECK(internal::get_scalar<std::int64_t>(table->GetColumnByName("bin2_id"), i) ==
        static_cast<std::int64_t>(p.bin2_id));

  CHECK(internal::get_scalar<N>(table->GetColumnByName("count"), i) == p.count);
}

template <std::int64_t i, typename N>
static void compare_pixel(const std::shared_ptr<arrow::Table>& table, const Pixel<N>& p) {
  assert(!!table);

  REQUIRE(i < table->num_rows());

  CHECK(internal::get_scalar<std::string>(table->GetColumnByName("chrom1"), i) ==
        p.coords.bin1.chrom().name());
  CHECK(internal::get_scalar<std::int32_t>(table->GetColumnByName("start1"), i) ==
        static_cast<std::int32_t>(p.coords.bin1.start()));
  CHECK(internal::get_scalar<std::int32_t>(table->GetColumnByName("end1"), i) ==
        static_cast<std::int32_t>(p.coords.bin1.end()));
  CHECK(internal::get_scalar<std::string>(table->GetColumnByName("chrom2"), i) ==
        p.coords.bin2.chrom().name());
  CHECK(internal::get_scalar<std::int32_t>(table->GetColumnByName("start2"), i) ==
        static_cast<std::int32_t>(p.coords.bin2.start()));
  CHECK(internal::get_scalar<std::int32_t>(table->GetColumnByName("end2"), i) ==
        static_cast<std::int32_t>(p.coords.bin2.end()));
  CHECK(internal::get_scalar<N>(table->GetColumnByName("count"), i) == p.count);
}

namespace internal {

[[nodiscard]] static std::vector<ThinPixel<std::uint8_t>> arrow_table_to_coo_vector(
    const std::shared_ptr<arrow::Table>& data) {
  assert(!!data);

  std::vector<ThinPixel<std::uint8_t>> buff(static_cast<std::size_t>(data->num_rows()));

  const auto bin1_ids = data->GetColumnByName("bin1_id");
  const auto bin2_ids = data->GetColumnByName("bin2_id");

  for (std::int64_t i = 0; i < data->num_rows(); ++i) {
    buff[static_cast<std::size_t>(i)].bin1_id =
        static_cast<std::uint64_t>(get_scalar<std::int64_t>(bin1_ids, i));
    buff[static_cast<std::size_t>(i)].bin2_id =
        static_cast<std::uint64_t>(get_scalar<std::int64_t>(bin2_ids, i));
  }
  return buff;
}

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
    buff[static_cast<std::size_t>(i)] =
        Pixel{chroms.at(get_scalar<std::string>(chrom1_ids, i)),
              static_cast<std::uint32_t>(get_scalar<std::int32_t>(start1, i)),
              static_cast<std::uint32_t>(get_scalar<std::int32_t>(end1, i)),
              chroms.at(get_scalar<std::string>(chrom2_ids, i)),
              static_cast<std::uint32_t>(get_scalar<std::int32_t>(start2, i)),
              static_cast<std::uint32_t>(get_scalar<std::int32_t>(end2, i)),
              std::uint8_t{}};
  }
  return buff;
}
}  // namespace internal

template <DataFrameFormat format, QuerySpan span>
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

template <DataFrameFormat format, QuerySpan span>
static void validate_diagonal_band(const Reference& chroms,
                                   const std::shared_ptr<arrow::Table>& table,
                                   std::uint64_t diagonal_band_width) {
  if constexpr (format == DataFrameFormat::COO) {
    const auto pixels = internal::arrow_table_to_coo_vector(table);
    if constexpr (span == QuerySpan::upper_triangle) {
      for (const auto& pixel : pixels) {
        CHECK(pixel.bin2_id - pixel.bin1_id < diagonal_band_width);
      }
    } else if constexpr (span == QuerySpan::lower_triangle) {
      for (const auto& pixel : pixels) {
        CHECK(pixel.bin1_id - pixel.bin2_id < diagonal_band_width);
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
      CHECK(pixel.coords.bin2.id() - pixel.coords.bin1.id() < diagonal_band_width);
    }
  } else if constexpr (span == QuerySpan::lower_triangle) {
    for (const auto& pixel : pixels) {
      CHECK(pixel.coords.bin1.id() - pixel.coords.bin2.id() < diagonal_band_width);
    }
  } else {
    throw std::logic_error("not implemented");
  }
}

TEST_CASE("Transformers (cooler): to dataframe", "[transformers][short]") {
  using N = std::int32_t;
  const auto path = datadir / "cooler" / "ENCFF993FGR.2500000.cool";
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

  SECTION("COO<int> upper_triangle (storage-mode=square)") {
    constexpr auto format = DataFrameFormat::COO;
    constexpr auto span = QuerySpan::upper_triangle;

    const cooler::File clr_square(
        (datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());
    const auto sel_ = clr_square.fetch();
    const auto table = ToDataFrame(sel_, sel_.begin<std::int32_t>(), format, nullptr, span)();

    CHECK(table->num_columns() == 3);
    CHECK(table->num_rows() == 53'154);
    CHECK(*table->column(2)->type() == *arrow::int32());

    // check head
    compare_pixel<0>(table, ThinPixel<N>{0, 0, 11768});
    compare_pixel<1>(table, ThinPixel<N>{0, 1, 14044});
    compare_pixel<2>(table, ThinPixel<N>{0, 2, 14496});

    // check tail
    compare_pixel<53151>(table, ThinPixel<N>{378, 378, 14432});
    compare_pixel<53152>(table, ThinPixel<N>{378, 379, 7150});
    compare_pixel<53153>(table, ThinPixel<N>{379, 379, 3534});

    validate_format<format, span>(clr.chromosomes(), table);
  }

  SECTION("COO<int> lower_triangle (storage-mode=square)") {
    constexpr auto format = DataFrameFormat::COO;
    constexpr auto span = QuerySpan::lower_triangle;

    const cooler::File clr_square(
        (datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());
    const auto sel_ = clr_square.fetch();
    const auto table =
        ToDataFrame(sel_, sel_.begin<std::int32_t>(), format, clr_square.bins_ptr(), span)();

    CHECK(table->num_columns() == 3);
    CHECK(table->num_rows() == 43'280);
    CHECK(*table->column(2)->type() == *arrow::int32());

    // check head
    compare_pixel<0>(table, ThinPixel<N>{0, 0, 11768});
    compare_pixel<1>(table, ThinPixel<N>{1, 0, 14081});
    compare_pixel<2>(table, ThinPixel<N>{1, 1, 14476});

    // check tail
    compare_pixel<43277>(table, ThinPixel<N>{379, 377, 6152});
    compare_pixel<43278>(table, ThinPixel<N>{379, 378, 7251});
    compare_pixel<43279>(table, ThinPixel<N>{379, 379, 3534});

    validate_format<format, span>(clr.chromosomes(), table);
  }

  SECTION("COO<int> full (storage-mode=square)") {
    constexpr auto format = DataFrameFormat::COO;
    constexpr auto span = QuerySpan::full;

    const cooler::File clr_square(
        (datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());
    const auto sel_ = clr_square.fetch();
    const auto table =
        ToDataFrame(sel_, sel_.begin<std::int32_t>(), format, clr_square.bins_ptr(), span)();

    CHECK(table->num_columns() == 3);
    CHECK(table->num_rows() == 96'133);
    CHECK(*table->column(2)->type() == *arrow::int32());

    // check head
    compare_pixel<0>(table, ThinPixel<N>{0, 0, 11768});
    compare_pixel<1>(table, ThinPixel<N>{0, 1, 14044});
    compare_pixel<2>(table, ThinPixel<N>{0, 2, 14496});

    // check tail
    compare_pixel<96130>(table, ThinPixel<N>{379, 377, 6152});
    compare_pixel<96131>(table, ThinPixel<N>{379, 378, 7251});
    compare_pixel<96132>(table, ThinPixel<N>{379, 379, 3534});
  }

  SECTION("BG2<int> upper_triangle (storage-mode=square)") {
    constexpr auto format = DataFrameFormat::BG2;
    constexpr auto span = QuerySpan::upper_triangle;

    const cooler::File clr_square(
        (datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());
    const auto& bins_square = clr_square.bins();
    const auto sel_ = clr_square.fetch();
    const auto table =
        ToDataFrame(sel_, sel_.begin<std::int32_t>(), format, clr_square.bins_ptr(), span)();

    CHECK(table->num_columns() == 7);
    CHECK(table->num_rows() == 53'154);
    CHECK(*table->column(6)->type() == *arrow::int32());

    // check head_square
    compare_pixel<0>(table, Pixel<N>{bins_square.at("chr1", 0), bins_square.at("chr1", 0), 11768});
    compare_pixel<1>(table,
                     Pixel<N>{bins_square.at("chr1", 0), bins_square.at("chr1", 8000), 14044});
    compare_pixel<2>(table,
                     Pixel<N>{bins_square.at("chr1", 0), bins_square.at("chr1", 16000), 14496});

    // check tail
    compare_pixel<53151>(
        table, Pixel<N>{bins_square.at("chr10", 288'000), bins_square.at("chr10", 288'000), 14432});
    compare_pixel<53152>(
        table, Pixel<N>{bins_square.at("chr10", 288'000), bins_square.at("chr10", 296'000), 7150});
    compare_pixel<53153>(
        table, Pixel<N>{bins_square.at("chr10", 296'000), bins_square.at("chr10", 296'000), 3534});

    validate_format<format, span>(clr.chromosomes(), table);
  }

  SECTION("BG2<int> lower_triangle (storage-mode=square)") {
    constexpr auto format = DataFrameFormat::BG2;
    constexpr auto span = QuerySpan::lower_triangle;

    const cooler::File clr_square(
        (datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());
    const auto& bins_square = clr_square.bins();
    const auto sel_ = clr_square.fetch();
    const auto table =
        ToDataFrame(sel_, sel_.begin<std::int32_t>(), format, clr_square.bins_ptr(), span)();

    CHECK(table->num_columns() == 7);
    CHECK(table->num_rows() == 43'280);
    CHECK(*table->column(6)->type() == *arrow::int32());

    // check head
    compare_pixel<0>(table, Pixel<N>{bins_square.at("chr1", 0), bins_square.at("chr1", 0), 11768});
    compare_pixel<1>(table,
                     Pixel<N>{bins_square.at("chr1", 8000), bins_square.at("chr1", 0), 14081});
    compare_pixel<2>(table,
                     Pixel<N>{bins_square.at("chr1", 8000), bins_square.at("chr1", 8000), 14476});

    // check tail
    compare_pixel<43277>(
        table, Pixel<N>{bins_square.at("chr10", 296'000), bins_square.at("chr10", 280'000), 6152});
    compare_pixel<43278>(
        table, Pixel<N>{bins_square.at("chr10", 296'000), bins_square.at("chr10", 288'000), 7251});
    compare_pixel<43279>(
        table, Pixel<N>{bins_square.at("chr10", 296'000), bins_square.at("chr10", 296'000), 3534});

    validate_format<format, span>(clr.chromosomes(), table);
  }

  SECTION("BG2<int> full (storage-mode=square)") {
    constexpr auto format = DataFrameFormat::BG2;
    constexpr auto span = QuerySpan::full;

    const cooler::File clr_square(
        (datadir / "cooler" / "cooler_storage_mode_square_test_file.mcool::/resolutions/8000")
            .string());
    const auto& bins_square = clr_square.bins();
    const auto sel_ = clr_square.fetch();
    const auto table =
        ToDataFrame(sel_, sel_.begin<std::int32_t>(), format, clr_square.bins_ptr(), span)();

    CHECK(table->num_columns() == 7);
    CHECK(table->num_rows() == 96'133);
    CHECK(*table->column(6)->type() == *arrow::int32());

    // check head
    compare_pixel<0>(table, Pixel<N>{bins_square.at("chr1", 0), bins_square.at("chr1", 0), 11768});
    compare_pixel<1>(table,
                     Pixel<N>{bins_square.at("chr1", 0), bins_square.at("chr1", 8000), 14044});
    compare_pixel<2>(table,
                     Pixel<N>{bins_square.at("chr1", 0), bins_square.at("chr1", 16000), 14496});

    // check tail
    compare_pixel<96130>(
        table, Pixel<N>{bins_square.at("chr10", 296'000), bins_square.at("chr10", 280'000), 6152});
    compare_pixel<96131>(
        table, Pixel<N>{bins_square.at("chr10", 296'000), bins_square.at("chr10", 288'000), 7251});
    compare_pixel<96132>(
        table, Pixel<N>{bins_square.at("chr10", 296'000), bins_square.at("chr10", 296'000), 3534});
  }

  SECTION("COO<int> upper_triangle w/ diagonal_band_width") {
    constexpr auto format = DataFrameFormat::COO;
    constexpr auto span = QuerySpan::upper_triangle;
    constexpr std::uint64_t diagonal_band_width{10};

    const auto table = ToDataFrame(first, last, format, nullptr, span, false, true, 256'000,
                                   diagonal_band_width)();

    CHECK(table->num_columns() == 3);
    CHECK(table->num_rows() == 856);
    CHECK(*table->column(2)->type() == *arrow::int32());

    // check head
    compare_pixel<0>(table, ThinPixel<std::int32_t>{0, 0, 266106});
    compare_pixel<1>(table, ThinPixel<std::int32_t>{0, 1, 32868});
    compare_pixel<2>(table, ThinPixel<std::int32_t>{0, 2, 13241});

    // check tail
    compare_pixel<853>(table, ThinPixel<std::int32_t>{98, 98, 1001844});
    compare_pixel<854>(table, ThinPixel<std::int32_t>{98, 99, 68621});
    compare_pixel<855>(table, ThinPixel<std::int32_t>{99, 99, 571144});

    validate_format<format, span>(clr.chromosomes(), table);
    validate_diagonal_band<format, span>(clr.chromosomes(), table, diagonal_band_width);
  }

  SECTION("BG2<int> upper_triangle w/ diagonal_band_width") {
    constexpr auto format = DataFrameFormat::BG2;
    constexpr auto span = QuerySpan::upper_triangle;
    constexpr std::uint64_t diagonal_band_width{10};

    const auto table = ToDataFrame(first, last, format, clr.bins_ptr(), span, false, true, 256'000,
                                   diagonal_band_width)();

    CHECK(table->num_columns() == 7);
    CHECK(table->num_rows() == 856);
    CHECK(*table->column(6)->type() == *arrow::int32());

    validate_format<format, span>(clr.chromosomes(), table);
    validate_diagonal_band<format, span>(clr.chromosomes(), table, diagonal_band_width);
  }

  SECTION("empty range") {
    const auto table = ToDataFrame(last, last)();
    CHECK(table->num_rows() == 0);
  }

  SECTION("invalid args") {
    const auto gw_sel = clr.fetch();

    constexpr auto format = DataFrameFormat::COO;
    constexpr auto span = QuerySpan::upper_triangle;
    constexpr std::uint64_t diagonal_band_width{10};

    const auto first_gw = gw_sel.begin<std::int32_t>();
    const auto last_gw = gw_sel.end<std::int32_t>();

    CHECK_THROWS_WITH(
        ToDataFrame(first_gw, last_gw, format, nullptr, span, false, true, 256'000,
                    diagonal_band_width)(),
        Catch::Matchers::ContainsSubstring("ToDataFrame<PixelIt>(): file index not loaded!"));
  }
}

TEST_CASE("Transformers (hic): to dataframe", "[transformers][short]") {
  const auto path = (datadir / "hic" / "4DNFIZ1ZVXC8.hic8").string();

  SECTION("ToDataFrame") {
    using N = std::int32_t;
    const hic::File hf(path, 2'500'000);
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

// NOLINTEND(*-avoid-magic-numbers, readability-function-cognitive-complexity)

}  // namespace hictk::test::transformers

#endif
