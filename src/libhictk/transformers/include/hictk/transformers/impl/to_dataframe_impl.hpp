// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <arrow/array.h>
#include <arrow/builder.h>
#include <arrow/table.h>
#include <arrow/type.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk::transformers {

template <typename PixelIt>
inline ToDataFrame<PixelIt>::ToDataFrame(PixelIt first, PixelIt last,
                                         std::shared_ptr<const BinTable> bins,
                                         std::size_t chunk_size)
    : _first(std::move(first)), _last(std::move(last)), _bins(std::move(bins)) {
  if (_bins) {
    _chrom_id_offset = _bins->chromosomes().at(0).is_all();
    const auto dict = make_chrom_dict(_bins->chromosomes());

    auto status = _chrom1_builder.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }

    status = _chrom2_builder.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
  }

  if (_bins) {
    _chrom1_id_buff.reserve(chunk_size);
    _start1_buff.reserve(chunk_size);
    _end1_buff.reserve(chunk_size);

    _chrom2_id_buff.reserve(chunk_size);
    _start2_buff.reserve(chunk_size);
    _end2_buff.reserve(chunk_size);
  } else {
    _bin1_id_buff.reserve(chunk_size);
    _bin2_id_buff.reserve(chunk_size);
  }

  _count_buff.reserve(chunk_size);
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::operator()() {
  if (_bins) {
    std::for_each(_first, _last, [&](const auto& p) { append(Pixel<N>(*_bins, p)); });
  } else {
    std::for_each(_first, _last, [&](const auto& p) { append(p); });
  }

  if (!_bin1_id.empty() || !_bin1_id_buff.empty()) {
    return make_coo_table();
  }
  return make_bg2_table();
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Schema> ToDataFrame<PixelIt>::coo_schema() const {
  return arrow::schema({
      // clang-format off
      arrow::field("bin1_id", arrow::uint64(), false),
      arrow::field("bin2_id", arrow::uint64(), false),
      arrow::field("count",   _count_builder.type(), false)
      // clang-format on
  });
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Schema> ToDataFrame<PixelIt>::bg2_schema() const {
  return arrow::schema({
      // clang-format off
      arrow::field("chrom1", arrow::dictionary(arrow::uint32(), arrow::utf8(), true), false),
      arrow::field("start1", arrow::uint32(), false),
      arrow::field("end1",   arrow::uint32(), false),
      arrow::field("chrom2", arrow::dictionary(arrow::uint32(), arrow::utf8(), true), false),
      arrow::field("start2", arrow::uint32(), false),
      arrow::field("end2",   arrow::uint32(), false),
      arrow::field("count",  _count_builder.type(), false)
      // clang-format on
  });
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append(const Pixel<N>& p) {
  if (_chrom1_id_buff.size() == _chrom1_id_buff.capacity()) {
    write_pixels();
  }

  _chrom1_id_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.chrom().id()) -
                            _chrom_id_offset);
  _start1_buff.push_back(p.coords.bin1.start());
  _end1_buff.push_back(p.coords.bin1.end());

  _chrom2_id_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.chrom().id()) -
                            _chrom_id_offset);
  _start2_buff.push_back(p.coords.bin2.start());
  _end2_buff.push_back(p.coords.bin2.end());

  _count_buff.push_back(p.count);
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append(const ThinPixel<N>& p) {
  if (_bin1_id_buff.size() == _bin1_id_buff.capacity()) {
    write_thin_pixels();
  }

  _bin1_id_buff.push_back(p.bin1_id);
  _bin2_id_buff.push_back(p.bin2_id);
  _count_buff.push_back(p.count);
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append(arrow::StringBuilder& builder, std::string_view data) {
  const auto status = builder.Append(std::string{data});
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
}

template <typename PixelIt>
template <typename ArrayBuilder, typename T>
inline void ToDataFrame<PixelIt>::append(ArrayBuilder& builder, const std::vector<T>& data) {
  const auto status = builder.AppendValues(data);
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append(arrow::StringDictionary32Builder& builder,
                                         const std::vector<std::int32_t>& data) {
  const auto status = builder.AppendIndices(data.data(), static_cast<std::int64_t>(data.size()));
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
}

template <typename PixelIt>
template <typename ArrayBuilder>
inline std::shared_ptr<arrow::Array> ToDataFrame<PixelIt>::finish(ArrayBuilder& builder) {
  auto result = builder.Finish();
  if (!result.status().ok()) {
    throw std::runtime_error(result.status().ToString());
  }

  return result.MoveValueUnsafe();
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::make_coo_table() {
  if (!_bin1_id_buff.empty()) {
    write_thin_pixels();
  }

  auto table = arrow::Table::Make(coo_schema(), {std::make_shared<arrow::ChunkedArray>(_bin1_id),
                                                 std::make_shared<arrow::ChunkedArray>(_bin2_id),
                                                 std::make_shared<arrow::ChunkedArray>(_count)});

  _bin1_id.clear();
  _bin2_id.clear();
  _count.clear();

  return table;
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::make_bg2_table() {
  if (!_chrom1_id_buff.empty()) {
    write_pixels();
  }

  // clang-format off
  auto table = arrow::Table::Make(
      bg2_schema(),
      {std::make_shared<arrow::ChunkedArray>(_chrom1),
       std::make_shared<arrow::ChunkedArray>(_start1),
       std::make_shared<arrow::ChunkedArray>(_end1),
       std::make_shared<arrow::ChunkedArray>(_chrom2),
       std::make_shared<arrow::ChunkedArray>(_start2),
       std::make_shared<arrow::ChunkedArray>(_end2),
       std::make_shared<arrow::ChunkedArray>(_count)});
  // clang-format on

  _chrom1.clear();
  _start1.clear();
  _end1.clear();

  _chrom2.clear();
  _start2.clear();
  _end2.clear();

  _count.clear();

  return table;
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Array> ToDataFrame<PixelIt>::make_chrom_dict(
    const hictk::Reference& chroms) {
  arrow::StringBuilder builder{};
  for (const auto& chrom : chroms) {
    if (!chrom.is_all()) {
      append(builder, std::string{chrom.name()});
    }
  }

  return finish(builder);
}

template <typename PixelIt>
void ToDataFrame<PixelIt>::write_thin_pixels() {
  if (!_bin1_id_buff.empty()) {
    append(_bin1_id_builder, _bin1_id_buff);
    append(_bin2_id_builder, _bin2_id_buff);
    append(_count_builder, _count_buff);

    _bin1_id.emplace_back(finish(_bin1_id_builder));
    _bin2_id.emplace_back(finish(_bin2_id_builder));
    _count.emplace_back(finish(_count_builder));

    _bin1_id_buff.clear();
    _bin2_id_buff.clear();
    _count_buff.clear();
  }
}

template <typename PixelIt>
void ToDataFrame<PixelIt>::write_pixels() {
  if (!_chrom1_id_buff.empty()) {
    append(_chrom1_builder, _chrom1_id_buff);
    append(_start1_builder, _start1_buff);
    append(_end1_builder, _end1_buff);

    append(_chrom2_builder, _chrom2_id_buff);
    append(_start2_builder, _start2_buff);
    append(_end2_builder, _end2_buff);

    append(_count_builder, _count_buff);

    _chrom1.emplace_back(finish(_chrom1_builder));
    _start1.emplace_back(finish(_start1_builder));
    _end1.emplace_back(finish(_end1_builder));

    _chrom2.emplace_back(finish(_chrom2_builder));
    _start2.emplace_back(finish(_start2_builder));
    _end2.emplace_back(finish(_end2_builder));

    _count.emplace_back(finish(_count_builder));

    _chrom1_id_buff.clear();
    _start1_buff.clear();
    _end1_buff.clear();

    _chrom2_id_buff.clear();
    _start2_buff.clear();
    _end2_buff.clear();

    _count_buff.clear();
  }
}

}  // namespace hictk::transformers
