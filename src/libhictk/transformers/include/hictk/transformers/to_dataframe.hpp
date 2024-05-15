// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_ARROW

#include <arrow/builder.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <arrow/type_traits.h>

#include <memory>
#include <type_traits>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"

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
  using N = decltype(std::declval<PixelIt>()->count);

  PixelIt _first{};
  PixelIt _last{};
  std::shared_ptr<const BinTable> _bins{};

  arrow::UInt64Builder _bin1_id{};
  arrow::UInt64Builder _bin2_id{};

  arrow::StringBuilder _chrom1{};
  arrow::UInt32Builder _start1{};
  arrow::UInt32Builder _end1{};

  arrow::StringBuilder _chrom2{};
  arrow::UInt32Builder _start2{};
  arrow::UInt32Builder _end2{};

  using NBuilder = decltype(internal::map_cpp_type_to_arrow_builder<N>());
  NBuilder _count{};

 public:
  ToDataFrame(PixelIt first_pixel, PixelIt last_pixel,
              std::shared_ptr<const BinTable> bins = nullptr);

  [[nodiscard]] std::shared_ptr<const arrow::Table> operator()();

 private:
  [[nodiscard]] std::shared_ptr<arrow::Schema> coo_schema() const;
  [[nodiscard]] std::shared_ptr<arrow::Schema> bg2_schema() const;

  void append(const Pixel<N>& p);
  void append(const ThinPixel<N>& p);

  template <typename ArrayBuilder, typename T>
  void append(ArrayBuilder& builder, const T& data);

  template <typename ArrayBuilder>
  [[nodiscard]] std::shared_ptr<arrow::Array> finish(ArrayBuilder& builder);

  [[nodiscard]] std::shared_ptr<const arrow::Table> make_coo_table();
  [[nodiscard]] std::shared_ptr<const arrow::Table> make_bg2_table();
};

}  // namespace hictk::transformers

#include "./impl/to_dataframe_impl.hpp"

#endif
