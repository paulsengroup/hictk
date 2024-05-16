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
                                         std::shared_ptr<const BinTable> bins)
    : _first(std::move(first)), _last(std::move(last)), _bins(std::move(bins)) {
  if (_bins) {
    const auto dict = make_chrom_dict(_bins->chromosomes());

    auto status = _chrom1.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }

    status = _chrom2.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }
  }
}
template <typename PixelIt>
inline std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::operator()() {
  if (_bins) {
    std::for_each(_first, _last, [&](const auto& p) { append(Pixel<N>(*_bins, p)); });
  } else {
    std::for_each(_first, _last, [&](const auto& p) { append(p); });
  }

  if (_bin1_id.length() != 0) {
    return make_coo_table();
  }
  return make_bg2_table();
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Schema> ToDataFrame<PixelIt>::coo_schema() const {
  return arrow::schema({
      // clang-format off
      arrow::field("bin1_id", arrow::uint64()),
      arrow::field("bin2_id", arrow::uint64()),
      arrow::field("count",   _count.type())
      // clang-format on
  });
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Schema> ToDataFrame<PixelIt>::bg2_schema() const {
  return arrow::schema({
      // clang-format off
      arrow::field("chrom1", arrow::dictionary(arrow::uint32(), arrow::utf8(), true)),
      arrow::field("start1", arrow::uint32()),
      arrow::field("end1",   arrow::uint32()),
      arrow::field("chrom2", arrow::dictionary(arrow::uint32(), arrow::utf8(), true)),
      arrow::field("start2", arrow::uint32()),
      arrow::field("end2",   arrow::uint32()),
      arrow::field("count",  _count.type())
      // clang-format on
  });
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append(const Pixel<N>& p) {
  append(_chrom1, std::string{p.coords.bin1.chrom().name()});
  append(_start1, p.coords.bin1.start());
  append(_end1, p.coords.bin1.end());

  append(_chrom2, std::string{p.coords.bin2.chrom().name()});
  append(_start2, p.coords.bin2.start());
  append(_end2, p.coords.bin2.end());

  append(_count, p.count);
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append(const ThinPixel<N>& p) {
  append(_bin1_id, p.bin1_id);
  append(_bin2_id, p.bin2_id);

  append(_count, p.count);
}

template <typename PixelIt>
template <typename ArrayBuilder, typename T>
inline void ToDataFrame<PixelIt>::append(ArrayBuilder& builder, const T& data) {
  const auto status = builder.Append(data);
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
  return arrow::Table::Make(coo_schema(), {finish(_bin1_id), finish(_bin2_id), finish(_count)});
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::make_bg2_table() {
  return arrow::Table::Make(bg2_schema(),
                            {finish(_chrom1), finish(_start1), finish(_end1), finish(_chrom2),
                             finish(_start2), finish(_end2), finish(_count)});
}

template <typename PixelIt>
std::shared_ptr<arrow::Array> ToDataFrame<PixelIt>::make_chrom_dict(
    const hictk::Reference& chroms) {
  arrow::StringBuilder builder{};
  for (const auto& chrom : chroms) {
    if (!chrom.is_all()) {
      append(builder, std::string{chrom.name()});
    }
  }

  return finish(builder);
}

}  // namespace hictk::transformers
