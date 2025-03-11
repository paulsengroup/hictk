// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// clang-format: off
#include "hictk/suppress_warnings.hpp"

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <arrow/array.h>
#include <arrow/builder.h>
#include <arrow/compute/api_vector.h>
#include <arrow/datum.h>
#include <arrow/table.h>
#include <arrow/type.h>
HICTK_DISABLE_WARNING_POP
// clang-format: on

#include <algorithm>
#include <cassert>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/transformers/common.hpp"
#include "hictk/transformers/diagonal_band.hpp"

namespace hictk::transformers {

template <typename PixelIt>
inline ToDataFrame<PixelIt>::ToDataFrame(PixelIt first, PixelIt last, DataFrameFormat format,
                                         std::shared_ptr<const BinTable> bins, QuerySpan span,
                                         bool include_bin_ids, bool mirror_pixels,
                                         std::size_t chunk_size,
                                         std::optional<std::uint64_t> diagonal_band_width)
    : _first(std::move(first)),
      _last(std::move(last)),
      _bins(std::move(bins)),
      _format(format),
      _span(span),
      _drop_bin_ids(!include_bin_ids),
      _mirror_pixels(mirror_pixels),
      _diagonal_band_width(diagonal_band_width),
      _chunk_size(chunk_size) {
  assert(chunk_size != 0);

  if (_format == DataFrameFormat::BG2 && !_bins) {
    throw std::runtime_error(
        "hictk::transformers::ToDataFrame: a bin table is required when format is "
        "DataFrameFormat::BG2");
  }

  if (_span != QuerySpan::upper_triangle && !_bins) {
    throw std::runtime_error(
        "hictk::transformers::ToDataFrame: a bin table is required when span is not "
        "QuerySpan::upper_triangle");
  }

  if (_format == DataFrameFormat::BG2) {
    assert(!!_bins);
    _chrom_id_offset = static_cast<std::int32_t>(_bins->chromosomes().at(0).is_all());
    const auto dict = make_chrom_dict(_bins->chromosomes());

    auto status = _chrom1_builder.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }

    status = _chrom2_builder.InsertMemoValues(*dict);
    if (!status.ok()) {
      throw std::runtime_error(status.ToString());
    }

    _chrom1_id_buff.reserve(_chunk_size);
    _start1_buff.reserve(_chunk_size);
    _end1_buff.reserve(_chunk_size);

    _chrom2_id_buff.reserve(_chunk_size);
    _start2_buff.reserve(_chunk_size);
    _end2_buff.reserve(_chunk_size);
  }

  if (_format == DataFrameFormat::BG2 || _span != QuerySpan::upper_triangle) {
    assert(!!_bins);
    _bin1_id_buff.reserve(_chunk_size);
    _bin2_id_buff.reserve(_chunk_size);
  }

  _count_buff.reserve(_chunk_size);
}

template <typename PixelIt>
template <typename PixelSelector>  // NOLINTNEXTLINE(*-unnecessary-value-param)
inline ToDataFrame<PixelIt>::ToDataFrame(const PixelSelector& sel, [[maybe_unused]] PixelIt it,
                                         DataFrameFormat format,
                                         std::shared_ptr<const BinTable> bins, QuerySpan span,
                                         bool include_bin_ids, std::size_t chunk_size,
                                         std::optional<std::uint64_t> diagonal_band_width)
    : ToDataFrame(sel.template begin<N>(), sel.template end<N>(), format, std::move(bins), span,
                  include_bin_ids, pixel_selector_is_symmetric_upper(sel), chunk_size,
                  diagonal_band_width) {}

template <typename PixelIt>
inline std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::operator()() {
  if (_diagonal_band_width.has_value()) {
    try {
      const DiagonalBand band_sel{_first, _last, *_diagonal_band_width};
      read_pixels(band_sel.begin(), band_sel.end());
    } catch (const std::runtime_error& e) {
      constexpr std::string_view prefix{"DiagonalBand<PixelIt>(): "};
      std::string_view msg{e.what()};
      if (msg.find(prefix) != 0) {
        throw;
      }

      msg.remove_prefix(prefix.size());
      throw std::runtime_error(fmt::format(
          FMT_STRING("ToDataFrame<PixelIt>(): {}. This only applies when diagonal_band_width is "
                     "specified when constructing a ToDataFrame instance."),
          msg));
    }
  } else {
    read_pixels(_first, _last);
  }

  if (_format == DataFrameFormat::COO) {
    return make_coo_table();
  }
  return make_bg2_table();
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Schema> ToDataFrame<PixelIt>::coo_schema() const {
  return arrow::schema({
      // clang-format off
      arrow::field("bin1_id", arrow::int64(), false),
      arrow::field("bin2_id", arrow::int64(), false),
      arrow::field("count",   _count_builder.type(), false)
      // clang-format on
  });
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Schema> ToDataFrame<PixelIt>::bg2_schema(bool with_bin_ids) const {
  arrow::FieldVector fields{};
  if (with_bin_ids) {
    fields.emplace_back(arrow::field("bin1_id", arrow::int64(), false));
    fields.emplace_back(arrow::field("bin2_id", arrow::int64(), false));
  }

  auto chrom_dict = dictionary(arrow::int32(), arrow::utf8());

  fields.emplace_back(arrow::field("chrom1", chrom_dict, false));
  fields.emplace_back(arrow::field("start1", arrow::int32(), false));
  fields.emplace_back(arrow::field("end1", arrow::int32(), false));
  fields.emplace_back(arrow::field("chrom2", chrom_dict, false));
  fields.emplace_back(arrow::field("start2", arrow::int32(), false));
  fields.emplace_back(arrow::field("end2", arrow::int32(), false));
  fields.emplace_back(arrow::field("count", _count_builder.type(), false));

  return arrow::schema(fields);
}

template <typename PixelIt>
template <typename PixelIt_>
inline void ToDataFrame<PixelIt>::read_pixels(PixelIt_ first, PixelIt_ last) {
  if (HICTK_LIKELY(_mirror_pixels)) {
    if (_format == DataFrameFormat::BG2) {
      assert(_bins);
      std::for_each(std::move(first), std::move(last),
                    [&](const auto& p) { append_symmetric(Pixel<N>(*_bins, p)); });
    } else {
      std::for_each(std::move(first), std::move(last),
                    [&](auto p) { append_symmetric(std::move(p)); });
    }
    return;
  }

  if (_format == DataFrameFormat::BG2) {
    assert(_bins);
    std::for_each(std::move(first), std::move(last),
                  [&](const auto& p) { append_asymmetric(Pixel<N>(*_bins, p)); });
  } else {
    std::for_each(std::move(first), std::move(last),
                  [&](auto p) { append_asymmetric(std::move(p)); });
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append_symmetric(Pixel<N> p) {
  assert(_mirror_pixels);

  if (_chrom1_id_buff.size() >= _chunk_size) {
    write_pixels();
  }

  if (_span == QuerySpan::lower_triangle) {
    std::swap(p.coords.bin1, p.coords.bin2);
  }

  const auto chrom1_id = static_cast<std::int32_t>(p.coords.bin1.chrom().id()) - _chrom_id_offset;
  const auto chrom2_id = static_cast<std::int32_t>(p.coords.bin2.chrom().id()) - _chrom_id_offset;

  assert(chrom1_id >= 0);
  assert(chrom2_id >= 0);

  if (_format == DataFrameFormat::BG2 || _span != QuerySpan::upper_triangle) {
    _bin1_id_buff.push_back(static_cast<std::int64_t>(p.coords.bin1.id()));
    _bin2_id_buff.push_back(static_cast<std::int64_t>(p.coords.bin2.id()));
  }

  _chrom1_id_buff.push_back(chrom1_id);
  _start1_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.start()));
  _end1_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.end()));

  _chrom2_id_buff.push_back(chrom2_id);
  _start2_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.start()));
  _end2_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.end()));

  _count_buff.push_back(p.count);

  if (_span == QuerySpan::full && p.coords.bin1 != p.coords.bin2) {
    _bin1_id_buff.push_back(static_cast<std::int64_t>(p.coords.bin2.id()));
    _bin2_id_buff.push_back(static_cast<std::int64_t>(p.coords.bin1.id()));

    _chrom1_id_buff.push_back(chrom2_id);
    _start1_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.start()));
    _end1_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.end()));

    _chrom2_id_buff.push_back(chrom1_id);
    _start2_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.start()));
    _end2_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.end()));

    _count_buff.push_back(p.count);
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append_symmetric(ThinPixel<N> p) {
  assert(_mirror_pixels);

  if (_bin1_id_buff.size() >= _chunk_size) {
    write_thin_pixels();
  }

  if (_span == QuerySpan::lower_triangle) {
    std::swap(p.bin1_id, p.bin2_id);
  }

  _bin1_id_buff.push_back(static_cast<std::int64_t>(p.bin1_id));
  _bin2_id_buff.push_back(static_cast<std::int64_t>(p.bin2_id));
  _count_buff.push_back(p.count);

  if (_span == QuerySpan::full && p.bin1_id != p.bin2_id) {
    _bin1_id_buff.push_back(static_cast<std::int64_t>(p.bin2_id));
    _bin2_id_buff.push_back(static_cast<std::int64_t>(p.bin1_id));
    _count_buff.push_back(p.count);
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append_asymmetric(const Pixel<N>& p) {
  assert(!_mirror_pixels);

  if (_chrom1_id_buff.size() >= _chunk_size) {
    write_pixels();
  }

  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  const auto populate_upper_triangle =
      _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  const auto chrom1_id = static_cast<std::int32_t>(p.coords.bin1.chrom().id()) - _chrom_id_offset;
  const auto chrom2_id = static_cast<std::int32_t>(p.coords.bin2.chrom().id()) - _chrom_id_offset;

  assert(chrom1_id >= 0);
  assert(chrom2_id >= 0);

  if (populate_upper_triangle && p.coords.bin1 <= p.coords.bin2) {
    if (_format == DataFrameFormat::BG2 || _span != QuerySpan::upper_triangle) {
      _bin1_id_buff.push_back(static_cast<std::int64_t>(p.coords.bin1.id()));
      _bin2_id_buff.push_back(static_cast<std::int64_t>(p.coords.bin2.id()));
    }

    _chrom1_id_buff.push_back(chrom1_id);
    _start1_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.start()));
    _end1_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.end()));

    _chrom2_id_buff.push_back(chrom2_id);
    _start2_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.start()));
    _end2_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.end()));

    _count_buff.push_back(p.count);
    return;
  }

  if (populate_lower_triangle && p.coords.bin1 >= p.coords.bin2) {
    if (_format == DataFrameFormat::BG2 || _span != QuerySpan::upper_triangle) {
      _bin1_id_buff.push_back(static_cast<std::int64_t>(p.coords.bin1.id()));
      _bin2_id_buff.push_back(static_cast<std::int64_t>(p.coords.bin2.id()));
    }

    _chrom1_id_buff.push_back(chrom1_id);
    _start1_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.start()));
    _end1_buff.push_back(static_cast<std::int32_t>(p.coords.bin1.end()));

    _chrom2_id_buff.push_back(chrom2_id);
    _start2_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.start()));
    _end2_buff.push_back(static_cast<std::int32_t>(p.coords.bin2.end()));

    _count_buff.push_back(p.count);
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append_asymmetric(ThinPixel<N> p) {
  assert(!_mirror_pixels);

  if (_bin1_id_buff.size() >= _chunk_size) {
    write_thin_pixels();
  }

  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  const auto populate_upper_triangle =
      _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  bool inserted = false;
  if (populate_upper_triangle && p.bin1_id <= p.bin2_id) {
    _bin1_id_buff.push_back(static_cast<std::int64_t>(p.bin1_id));
    _bin2_id_buff.push_back(static_cast<std::int64_t>(p.bin2_id));
    _count_buff.push_back(p.count);
    inserted = true;
  }

  if (populate_lower_triangle && !inserted && p.bin1_id >= p.bin2_id) {
    _bin1_id_buff.push_back(static_cast<std::int64_t>(p.bin1_id));
    _bin2_id_buff.push_back(static_cast<std::int64_t>(p.bin2_id));
    _count_buff.push_back(p.count);
  }
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

  if (_bin1_id.empty()) {
    assert(_bin2_id.empty());
    auto result = arrow::Table::MakeEmpty(coo_schema());
    if (!result.ok()) {
      throw std::runtime_error(result.status().ToString());
    }
    return result.MoveValueUnsafe();
  }

  auto table = arrow::Table::Make(coo_schema(), {std::make_shared<arrow::ChunkedArray>(_bin1_id),
                                                 std::make_shared<arrow::ChunkedArray>(_bin2_id),
                                                 std::make_shared<arrow::ChunkedArray>(_count)});
  assert(table->ValidateFull().ok());

  _bin1_id.clear();
  _bin2_id.clear();
  _count.clear();

  if (_span != QuerySpan::upper_triangle) {
    table = sort_table(table);
  }

  return table;
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::make_bg2_table() {
  if (!_chrom1_id_buff.empty()) {
    write_pixels();
  }

  if (_chrom1.empty()) {
    auto result = arrow::Table::MakeEmpty(bg2_schema(!_drop_bin_ids));
    if (!result.ok()) {
      throw std::runtime_error(result.status().ToString());
    }
    return result.MoveValueUnsafe();
  }

  std::shared_ptr<arrow::Table> table{};
  if (!_bin1_id.empty()) {
    assert(!_bin2_id.empty());

    // clang-format off
    table = arrow::Table::Make(
      bg2_schema(true),
      {std::make_shared<arrow::ChunkedArray>(_bin1_id),
       std::make_shared<arrow::ChunkedArray>(_bin2_id),
       std::make_shared<arrow::ChunkedArray>(_chrom1),
       std::make_shared<arrow::ChunkedArray>(_start1),
       std::make_shared<arrow::ChunkedArray>(_end1),
       std::make_shared<arrow::ChunkedArray>(_chrom2),
       std::make_shared<arrow::ChunkedArray>(_start2),
       std::make_shared<arrow::ChunkedArray>(_end2),
       std::make_shared<arrow::ChunkedArray>(_count)});
    // clang-format on
  } else {
    assert(_bin1_id.empty());
    assert(_bin2_id.empty());
    // clang-format off
    table = arrow::Table::Make(
      bg2_schema(),
      {std::make_shared<arrow::ChunkedArray>(_chrom1),
       std::make_shared<arrow::ChunkedArray>(_start1),
       std::make_shared<arrow::ChunkedArray>(_end1),
       std::make_shared<arrow::ChunkedArray>(_chrom2),
       std::make_shared<arrow::ChunkedArray>(_start2),
       std::make_shared<arrow::ChunkedArray>(_end2),
       std::make_shared<arrow::ChunkedArray>(_count)});
    // clang-format on
  }
  assert(table->ValidateFull().ok());

  _bin1_id.clear();
  _bin2_id.clear();

  _chrom1.clear();
  _start1.clear();
  _end1.clear();

  _chrom2.clear();
  _start2.clear();
  _end2.clear();

  _count.clear();

  if (_span != QuerySpan::upper_triangle) {
    table = sort_table(table);

    if (_drop_bin_ids) {
      auto result = table->RemoveColumn(0);
      if (!result.ok()) {
        throw std::runtime_error(result.status().ToString());
      }
      table = result.MoveValueUnsafe();
      result = table->RemoveColumn(0);
      if (!result.ok()) {
        throw std::runtime_error(result.status().ToString());
      }
      table = result.MoveValueUnsafe();
      assert(table->ValidateFull().ok());
    }
  }

  return table;
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Array> ToDataFrame<PixelIt>::make_chrom_dict(
    const Reference& chroms) {
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
  const auto populate_lower_triangle = _span != QuerySpan::upper_triangle;
  if (!_chrom1_id_buff.empty()) {
    if (populate_lower_triangle) {
      append(_bin1_id_builder, _bin1_id_buff);
      append(_bin2_id_builder, _bin2_id_buff);
    }

    append(_chrom1_builder, _chrom1_id_buff);
    append(_start1_builder, _start1_buff);
    append(_end1_builder, _end1_buff);

    append(_chrom2_builder, _chrom2_id_buff);
    append(_start2_builder, _start2_buff);
    append(_end2_builder, _end2_buff);

    append(_count_builder, _count_buff);

    if (populate_lower_triangle) {
      _bin1_id.emplace_back(finish(_bin1_id_builder));
      _bin2_id.emplace_back(finish(_bin2_id_builder));
    }

    _chrom1.emplace_back(finish(_chrom1_builder));
    _start1.emplace_back(finish(_start1_builder));
    _end1.emplace_back(finish(_end1_builder));

    _chrom2.emplace_back(finish(_chrom2_builder));
    _start2.emplace_back(finish(_start2_builder));
    _end2.emplace_back(finish(_end2_builder));

    _count.emplace_back(finish(_count_builder));

    _bin1_id_buff.clear();
    _bin2_id_buff.clear();

    _chrom1_id_buff.clear();
    _start1_buff.clear();
    _end1_buff.clear();

    _chrom2_id_buff.clear();
    _start2_buff.clear();
    _end2_buff.clear();

    _count_buff.clear();
  }
}

template <typename PixelIt>
std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::sort_table(
    std::shared_ptr<arrow::Table> table) {
  const arrow::compute::SortOptions opts{
      {arrow::compute::SortKey{"bin1_id", arrow::compute::SortOrder::Ascending},
       arrow::compute::SortKey{"bin2_id", arrow::compute::SortOrder::Ascending}}};

  const arrow::Datum vtable{std::move(table)};

  const auto arg_sorter = arrow::compute::SortIndices(vtable, opts);
  if (!arg_sorter.ok()) {
    throw std::runtime_error(arg_sorter.status().ToString());
  }

  if constexpr (ndebug_not_defined()) {
    return arrow::compute::Take(vtable, *arg_sorter)->table();
  } else {
    table = arrow::compute::Take(vtable, *arg_sorter, arrow::compute::TakeOptions::NoBoundsCheck())
                ->table();
    assert(table->ValidateFull().ok());
    return table;
  }
}

template <typename PixelIt>
template <typename PixelSelector>
inline bool ToDataFrame<PixelIt>::pixel_selector_is_symmetric_upper(
    const PixelSelector& sel) noexcept {
  if constexpr (std::is_same_v<PixelSelector, cooler::PixelSelector>) {
    return sel.is_symmetric_upper();
  }
  return true;
}

}  // namespace hictk::transformers
