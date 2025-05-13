// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_ARROW

// clang-format: off
#include "hictk/suppress_warnings.hpp"

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/builder.h>
#include <arrow/table.h>
#include <arrow/type.h>
HICTK_DISABLE_WARNING_POP
// clang-format: on

#include <cstdint>
#include <memory>
#include <optional>
#include <type_traits>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/transformers/common.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::transformers {

namespace internal {
template <typename N, typename std::enable_if_t<std::is_same_v<N, std::uint8_t>>* = nullptr>
static arrow::UInt8Builder map_cpp_type_to_arrow_builder() {
  return arrow::UInt8Builder{};
}
template <typename N, typename std::enable_if_t<std::is_same_v<N, std::uint16_t>>* = nullptr>
static arrow::UInt16Builder map_cpp_type_to_arrow_builder() {
  return arrow::UInt16Builder{};
}
template <typename N, typename std::enable_if_t<std::is_same_v<N, std::uint32_t>>* = nullptr>
static arrow::UInt32Builder map_cpp_type_to_arrow_builder() {
  return arrow::UInt32Builder{};
}
template <typename N, typename std::enable_if_t<std::is_same_v<N, std::uint64_t>>* = nullptr>
static arrow::UInt64Builder map_cpp_type_to_arrow_builder() {
  return arrow::UInt64Builder{};
}

template <typename N, typename std::enable_if_t<std::is_same_v<N, std::int8_t>>* = nullptr>
static arrow::Int8Builder map_cpp_type_to_arrow_builder() {
  return arrow::Int8Builder{};
}
template <typename N, typename std::enable_if_t<std::is_same_v<N, std::int16_t>>* = nullptr>
static arrow::Int16Builder map_cpp_type_to_arrow_builder() {
  return arrow::Int16Builder{};
}
template <typename N, typename std::enable_if_t<std::is_same_v<N, std::int32_t>>* = nullptr>
static arrow::Int32Builder map_cpp_type_to_arrow_builder() {
  return arrow::Int32Builder{};
}
template <typename N, typename std::enable_if_t<std::is_same_v<N, std::int64_t>>* = nullptr>
static arrow::Int64Builder map_cpp_type_to_arrow_builder() {
  return arrow::Int64Builder{};
}

template <typename N, typename std::enable_if_t<std::is_same_v<N, float>>* = nullptr>
static arrow::FloatBuilder map_cpp_type_to_arrow_builder() {
  return arrow::FloatBuilder{};
}
template <typename N, typename std::enable_if_t<std::is_same_v<N, double>>* = nullptr>
static arrow::DoubleBuilder map_cpp_type_to_arrow_builder() {
  return arrow::DoubleBuilder{};
}
}  // namespace internal

enum class DataFrameFormat : std::uint_fast8_t { COO, BG2 };

template <typename PixelIt>
class ToDataFrame {
  using PixelT = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  using N = decltype(std::declval<PixelIt>()->count);
  static_assert(std::is_same_v<PixelT, ThinPixel<N>>);

  struct Builder {
    arrow::StringDictionary32Builder chrom{};
    arrow::Int64Builder int64{};
    arrow::Int32Builder int32{};
    using NBuilder = decltype(internal::map_cpp_type_to_arrow_builder<N>());
    NBuilder count{};

    std::size_t chunk_size{};
    std::int32_t chrom_id_offset{};

    Builder() = default;
    Builder(const Reference& chroms, std::size_t chunk_size_);
    explicit Builder(std::size_t chunk_size_);

    [[nodiscard]] std::shared_ptr<arrow::DataType> count_type() const;

   private:
    [[nodiscard]] static std::shared_ptr<arrow::Array> make_chrom_dict(const Reference& chroms);
  };

  struct Buffer {
    std::vector<std::int64_t> bin1_id{};
    std::vector<std::int64_t> bin2_id{};

    std::vector<std::int32_t> chrom1_id{};
    std::vector<std::int32_t> start1{};
    std::vector<std::int32_t> end1{};

    std::vector<std::int32_t> chrom2_id{};
    std::vector<std::int32_t> start2{};
    std::vector<std::int32_t> end2{};

    std::vector<N> count{};

    Buffer() = default;
    Buffer(DataFrameFormat format, QuerySpan span, std::size_t chunk_size_);

    [[nodiscard]] std::size_t capacity() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool empty() const noexcept;

    void clear() noexcept;
  };

  struct VectorChunks {
    arrow::ArrayVector bin1_id{};
    arrow::ArrayVector bin2_id{};

    arrow::ArrayVector chrom1{};
    arrow::ArrayVector start1{};
    arrow::ArrayVector end1{};

    arrow::ArrayVector chrom2{};
    arrow::ArrayVector start2{};
    arrow::ArrayVector end2{};

    arrow::ArrayVector count{};

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool empty() const noexcept;

    void clear();
  };

  PixelIt _first{};
  PixelIt _last{};
  std::shared_ptr<const BinTable> _bins{};
  std::optional<PixelCoordinates> _coord1{};
  std::optional<PixelCoordinates> _coord2{};

  DataFrameFormat _format{DataFrameFormat::COO};
  QuerySpan _span{QuerySpan::upper_triangle};
  bool _drop_bin_ids{true};
  bool _mirror_pixels{true};
  std::optional<std::uint64_t> _diagonal_band_width{};

  std::unique_ptr<Builder> _builder{};
  std::unique_ptr<Buffer> _buffer{};
  std::unique_ptr<VectorChunks> _chunks{};

 public:
  // NOLINTBEGIN(*-avoid-magic-numbers, *-unnecessary-value-param)
  ToDataFrame(PixelIt first_pixel, PixelIt last_pixel, std::optional<PixelCoordinates> coord1_,
              std::optional<PixelCoordinates> coord2_,
              DataFrameFormat format = DataFrameFormat::COO,
              std::shared_ptr<const BinTable> bins = nullptr,
              QuerySpan span = QuerySpan::upper_triangle, bool include_bin_ids = false,
              bool mirror_pixels = true, std::size_t chunk_size = 256'000,
              std::optional<std::uint64_t> diagonal_band_width = {});

  ToDataFrame(PixelIt first_pixel, PixelIt last_pixel,
              DataFrameFormat format = DataFrameFormat::COO,
              std::shared_ptr<const BinTable> bins = nullptr,
              QuerySpan span = QuerySpan::upper_triangle, bool include_bin_ids = false,
              bool mirror_pixels = true, std::size_t chunk_size = 256'000,
              std::optional<std::uint64_t> diagonal_band_width = {});

  template <typename PixelSelector,
            typename std::enable_if_t<internal::has_coord1_member_fx<PixelSelector>>* = nullptr>
  ToDataFrame(const PixelSelector& sel, PixelIt it, DataFrameFormat format = DataFrameFormat::COO,
              std::shared_ptr<const BinTable> bins = nullptr,
              QuerySpan span = QuerySpan::upper_triangle, bool include_bin_ids = false,
              std::size_t chunk_size = 256'000,
              std::optional<std::uint64_t> diagonal_band_width = {});
  template <typename PixelSelector,
            typename std::enable_if_t<!internal::has_coord1_member_fx<PixelSelector>>* = nullptr>
  ToDataFrame(const PixelSelector& sel, PixelIt it, DataFrameFormat format = DataFrameFormat::COO,
              std::shared_ptr<const BinTable> bins = nullptr,
              QuerySpan span = QuerySpan::upper_triangle, bool include_bin_ids = false,
              std::size_t chunk_size = 256'000,
              std::optional<std::uint64_t> diagonal_band_width = {});
  // NOLINTEND(*-avoid-magic-numbers, *-unnecessary-value-param)

  [[nodiscard]] std::shared_ptr<arrow::Table> operator()();

 private:
  [[nodiscard]] std::shared_ptr<arrow::Schema> coo_schema() const;
  [[nodiscard]] std::shared_ptr<arrow::Schema> bg2_schema(bool with_bin_ids = false) const;

  template <typename PixelIt_>
  void read_pixels(PixelIt_ first, PixelIt_ last);

  void append_symmetric(Pixel<N> p);
  void append_symmetric(ThinPixel<N> p);

  void append_asymmetric(const Pixel<N>& p);
  void append_asymmetric(ThinPixel<N> p);

  void append(const Pixel<N>& p);
  void append(const ThinPixel<N>& p);

  [[nodiscard]] bool overlaps(const ThinPixel<N>& p) const noexcept;
  [[nodiscard]] bool overlaps(const Pixel<N>& p) const noexcept;
  [[nodiscard]] bool overlaps(std::uint64_t bin1_id, std::uint64_t bin2_id) const noexcept;

  static void append(arrow::StringBuilder& builder, std::string_view data);
  template <typename ArrayBuilder, typename T>
  static void append(ArrayBuilder& builder, const std::vector<T>& data);
  static void append(arrow::StringDictionary32Builder& builder,
                     const std::vector<std::int32_t>& data);

  template <typename ArrayBuilder>
  [[nodiscard]] static std::shared_ptr<arrow::Array> finish(ArrayBuilder& builder);

  [[nodiscard]] std::shared_ptr<arrow::Table> make_coo_table();
  [[nodiscard]] std::shared_ptr<arrow::Table> make_bg2_table();

  void commit_thin_pixels();
  void commit_pixels();

  template <typename PixelSelector>
  [[nodiscard]] static bool pixel_selector_is_symmetric_upper(const PixelSelector& sel) noexcept;

  static std::shared_ptr<arrow::Table> sort_table(std::shared_ptr<arrow::Table> table);

  [[nodiscard]] static QuerySpan fix_query_span(const std::optional<PixelCoordinates>& coord1,
                                                const std::optional<PixelCoordinates>& coord2,
                                                QuerySpan requested_span);
  [[nodiscard]] static std::optional<PixelCoordinates> fix_coordinates(
      std::optional<PixelCoordinates> coord);
};

}  // namespace hictk::transformers

#include "./impl/to_dataframe_impl.hpp"

#endif
