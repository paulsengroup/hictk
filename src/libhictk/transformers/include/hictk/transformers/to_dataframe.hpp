// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_ARROW

#include <arrow/array.h>
#include <arrow/builder.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <arrow/type_traits.h>

#include <memory>
#include <type_traits>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
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

template <typename PixelIt>
class ToDataFrame {
  using PixelT = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  using N = decltype(std::declval<PixelIt>()->count);
  static_assert(std::is_same_v<PixelT, hictk::ThinPixel<N>>);

  PixelIt _first{};
  PixelIt _last{};
  std::shared_ptr<const BinTable> _bins{};

  arrow::UInt64Builder _bin1_id_builder{};
  arrow::UInt64Builder _bin2_id_builder{};

  arrow::StringDictionary32Builder _chrom1_builder{};
  arrow::UInt32Builder _start1_builder{};
  arrow::UInt32Builder _end1_builder{};

  arrow::StringDictionary32Builder _chrom2_builder{};
  arrow::UInt32Builder _start2_builder{};
  arrow::UInt32Builder _end2_builder{};

  using NBuilder = decltype(internal::map_cpp_type_to_arrow_builder<N>());
  NBuilder _count_builder{};

  std::vector<std::uint64_t> _bin1_id_buff{};
  std::vector<std::uint64_t> _bin2_id_buff{};

  std::vector<std::int32_t> _chrom1_id_buff{};
  std::vector<std::uint32_t> _start1_buff{};
  std::vector<std::uint32_t> _end1_buff{};

  std::vector<std::int32_t> _chrom2_id_buff{};
  std::vector<std::uint32_t> _start2_buff{};
  std::vector<std::uint32_t> _end2_buff{};

  std::vector<N> _count_buff{};

  arrow::ArrayVector _bin1_id{};
  arrow::ArrayVector _bin2_id{};

  arrow::ArrayVector _chrom1{};
  arrow::ArrayVector _start1{};
  arrow::ArrayVector _end1{};

  arrow::ArrayVector _chrom2{};
  arrow::ArrayVector _start2{};
  arrow::ArrayVector _end2{};

  arrow::ArrayVector _count{};

  std::int32_t _chrom_id_offset{};

 public:
  ToDataFrame(PixelIt first_pixel, PixelIt last_pixel,
              std::shared_ptr<const BinTable> bins = nullptr, std::size_t chunk_size = 256'000);

  [[nodiscard]] std::shared_ptr<arrow::Table> operator()();

 private:
  [[nodiscard]] std::shared_ptr<arrow::Schema> coo_schema() const;
  [[nodiscard]] std::shared_ptr<arrow::Schema> bg2_schema() const;

  void append(const Pixel<N>& p);
  void append(const ThinPixel<N>& p);

  static void append(arrow::StringBuilder& builder, std::string_view data);
  template <typename ArrayBuilder, typename T>
  static void append(ArrayBuilder& builder, const std::vector<T>& data);
  static void append(arrow::StringDictionary32Builder& builder,
                     const std::vector<std::int32_t>& data);

  template <typename ArrayBuilder>
  [[nodiscard]] static std::shared_ptr<arrow::Array> finish(ArrayBuilder& builder);

  [[nodiscard]] std::shared_ptr<arrow::Table> make_coo_table();
  [[nodiscard]] std::shared_ptr<arrow::Table> make_bg2_table();

  [[nodiscard]] static std::shared_ptr<arrow::Array> make_chrom_dict(
      const hictk::Reference& chroms);

  void write_thin_pixels();
  void write_pixels();
};

}  // namespace hictk::transformers

#include "./impl/to_dataframe_impl.hpp"

#endif
