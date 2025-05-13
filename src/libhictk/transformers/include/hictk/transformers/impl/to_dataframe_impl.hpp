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
inline ToDataFrame<PixelIt>::Builder::Builder(std::size_t chunk_size_) : Builder({}, chunk_size_) {}

template <typename PixelIt>
inline ToDataFrame<PixelIt>::Builder::Builder(const Reference& chroms, std::size_t chunk_size_)
    : chunk_size(chunk_size_),
      chrom_id_offset(chroms.empty() ? std::int32_t{}
                                     : static_cast<std::int32_t>(chroms.at(0).is_all())) {
  if (chunk_size == 0) {
    throw std::invalid_argument("chunk size cannot be 0");
  }

  const auto dict = make_chrom_dict(chroms);
  if (!dict) {
    return;
  }

  const auto status = chrom.InsertMemoValues(*dict);
  if (!status.ok()) {
    throw std::runtime_error(status.ToString());
  }
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Array> ToDataFrame<PixelIt>::Builder::make_chrom_dict(
    const Reference& chroms) {
  if (chroms.remove_ALL().empty()) {
    return nullptr;
  }

  arrow::StringBuilder builder{};
  for (const auto& chrom : chroms) {
    if (!chrom.is_all()) {
      append(builder, std::string{chrom.name()});
    }
  }

  return finish(builder);
}

template <typename PixelIt>
inline std::shared_ptr<arrow::DataType> ToDataFrame<PixelIt>::Builder::count_type() const {
  return count.type();
}

template <typename PixelIt>
inline ToDataFrame<PixelIt>::Buffer::Buffer(DataFrameFormat format, QuerySpan span,
                                            std::size_t chunk_size_) {
  if (chunk_size_ == 0) {
    throw std::invalid_argument("chunk size cannot be 0");
  }

  if (format == DataFrameFormat::BG2) {
    chrom1_id.reserve(chunk_size_);
    start1.reserve(chunk_size_);
    end1.reserve(chunk_size_);

    chrom2_id.reserve(chunk_size_);
    start2.reserve(chunk_size_);
    end2.reserve(chunk_size_);
  }

  if (format == DataFrameFormat::COO || span != QuerySpan::upper_triangle) {
    bin1_id.reserve(chunk_size_);
    bin2_id.reserve(chunk_size_);
  }

  count.reserve(chunk_size_);
}

template <typename PixelIt>
inline std::size_t ToDataFrame<PixelIt>::Buffer::capacity() const noexcept {
  return count.capacity();
}

template <typename PixelIt>
inline std::size_t ToDataFrame<PixelIt>::Buffer::size() const noexcept {
  return count.size();
}

template <typename PixelIt>
inline bool ToDataFrame<PixelIt>::Buffer::empty() const noexcept {
  return size() == 0;
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::Buffer::clear() noexcept {
  bin1_id.clear();
  bin2_id.clear();

  chrom1_id.clear();
  start1.clear();
  end1.clear();

  chrom2_id.clear();
  start2.clear();
  end2.clear();

  count.clear();
}

template <typename PixelIt>
inline std::size_t ToDataFrame<PixelIt>::VectorChunks::size() const noexcept {
  return count.size();
}

template <typename PixelIt>
inline bool ToDataFrame<PixelIt>::VectorChunks::empty() const noexcept {
  return size() == 0;
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::VectorChunks::clear() {
  bin1_id.clear();
  bin2_id.clear();
  chrom1.clear();
  start1.clear();
  end1.clear();
  chrom2.clear();
  start2.clear();
  end2.clear();
  count.clear();
}

template <typename PixelIt>
inline ToDataFrame<PixelIt>::ToDataFrame(
    PixelIt first, PixelIt last, std::optional<PixelCoordinates> coord1_,
    std::optional<PixelCoordinates> coord2_, DataFrameFormat format,
    std::shared_ptr<const BinTable> bins, QuerySpan span, bool include_bin_ids, bool mirror_pixels,
    std::size_t chunk_size, std::optional<std::uint64_t> diagonal_band_width)
    : _first(std::move(first)),
      _last(std::move(last)),
      _bins(std::move(bins)),
      _coord1(fix_coordinates(std::move(coord1_))),
      _coord2(fix_coordinates(coord2_.has_value() ? std::move(coord2_) : _coord1)),
      _format(format),
      _span(mirror_pixels ? fix_query_span(_coord1, _coord2, span) : span),
      _drop_bin_ids(!include_bin_ids),
      _mirror_pixels(mirror_pixels),
      _diagonal_band_width(diagonal_band_width),
      _builder(std::make_unique<Builder>(_bins ? _bins->chromosomes() : Reference{}, chunk_size)),
      _buffer(std::make_unique<Buffer>(_format, _span, chunk_size)),
      _chunks(std::make_unique<VectorChunks>()) {
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
}

template <typename PixelIt>
inline ToDataFrame<PixelIt>::ToDataFrame(PixelIt first, PixelIt last, DataFrameFormat format,
                                         std::shared_ptr<const BinTable> bins, QuerySpan span,
                                         bool include_bin_ids, bool mirror_pixels,
                                         std::size_t chunk_size,
                                         std::optional<std::uint64_t> diagonal_band_width)
    : ToDataFrame(std::move(first), std::move(last), std::nullopt, std::nullopt, format,
                  std::move(bins), span, include_bin_ids, mirror_pixels, chunk_size,
                  diagonal_band_width) {}

template <typename PixelIt>
template <typename PixelSelector,  // NOLINTNEXTLINE(modernize-type-traits)
          typename std::enable_if_t<internal::has_coord1_member_fx<PixelSelector>>*>
inline ToDataFrame<PixelIt>::ToDataFrame(const PixelSelector& sel, [[maybe_unused]] PixelIt it,
                                         DataFrameFormat format,
                                         std::shared_ptr<const BinTable> bins, QuerySpan span,
                                         bool include_bin_ids, std::size_t chunk_size,
                                         std::optional<std::uint64_t> diagonal_band_width)
    : ToDataFrame(sel.template begin<N>(), sel.template end<N>(), sel.coord1(), sel.coord2(),
                  format, std::move(bins), span, include_bin_ids,
                  pixel_selector_is_symmetric_upper(sel), chunk_size, diagonal_band_width) {
  if constexpr (internal::has_coord1_member_fx<PixelSelector>) {
    if (_span != QuerySpan::full) {
      return;
    }

    if (!_coord1.has_value()) {
      assert(!_coord2.has_value());
      return;
    }

    // NOLINTBEGIN(*-unchecked-optional-access)
    if (_coord1->bin1.chrom() == _coord2->bin1.chrom() && *_coord1 != *_coord2) {
      auto bin1 = std::min(_coord1->bin1, _coord2->bin1);
      auto bin2 = std::max(_coord1->bin2, _coord2->bin2);

      const PixelCoordinates coord{std::move(bin1), std::move(bin2)};
      const auto new_sel = sel.fetch(coord, coord);
      *this =
          ToDataFrame(new_sel.template begin<N>(), new_sel.template end<N>(), std::move(_coord1),
                      std::move(_coord2), _format, std::move(_bins), _span, !_drop_bin_ids,
                      _mirror_pixels, _buffer->capacity(), _diagonal_band_width);
    }
    // NOLINTEND(*-unchecked-optional-access)
  }
}

template <typename PixelIt>
template <typename PixelSelector,  // NOLINTNEXTLINE(modernize-type-traits)
          typename std::enable_if_t<!internal::has_coord1_member_fx<PixelSelector>>*>
inline ToDataFrame<PixelIt>::ToDataFrame(const PixelSelector& sel, [[maybe_unused]] PixelIt it,
                                         DataFrameFormat format,
                                         std::shared_ptr<const BinTable> bins, QuerySpan span,
                                         bool include_bin_ids, std::size_t chunk_size,
                                         std::optional<std::uint64_t> diagonal_band_width)
    : ToDataFrame(sel.template begin<N>(), sel.template end<N>(), std::nullopt, std::nullopt,
                  format, std::move(bins), span, include_bin_ids,
                  pixel_selector_is_symmetric_upper(sel), chunk_size, diagonal_band_width) {}

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
  assert(_builder);
  return arrow::schema({
      // clang-format off
      arrow::field("bin1_id", arrow::int64(), false),
      arrow::field("bin2_id", arrow::int64(), false),
      arrow::field("count",   _builder->count_type(), false)
      // clang-format on
  });
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Schema> ToDataFrame<PixelIt>::bg2_schema(bool with_bin_ids) const {
  assert(_builder);

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
  fields.emplace_back(arrow::field("count", _builder->count_type(), false));

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

  auto do_append = [&](bool swap) {
    if (swap) {
      std::swap(p.coords.bin1, p.coords.bin2);
    }
    if (HICTK_LIKELY(overlaps(p))) {
      append(p);
    }
  };

  switch (_span) {
    case QuerySpan::upper_triangle: {
      do_append(false);
      break;
    }
    case QuerySpan::lower_triangle: {
      do_append(true);
      break;
    }
    case QuerySpan::full: {
      do_append(false);
      if (p.coords.bin1 != p.coords.bin2) {
        do_append(true);
      }
    }
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append_symmetric(ThinPixel<N> p) {
  assert(_mirror_pixels);
  assert(_mirror_pixels);

  auto do_append = [&](bool swap) {
    if (swap) {
      std::swap(p.bin1_id, p.bin2_id);
    }
    if (HICTK_LIKELY(overlaps(p))) {
      append(p);
    }
  };

  switch (_span) {
    case QuerySpan::upper_triangle: {
      do_append(false);
      break;
    }
    case QuerySpan::lower_triangle: {
      do_append(true);
      break;
    }
    case QuerySpan::full: {
      do_append(false);
      if (p.bin1_id != p.bin2_id) {
        do_append(true);
      }
    }
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append_asymmetric(const Pixel<N>& p) {
  assert(!_mirror_pixels);

  if (HICTK_UNLIKELY(!overlaps(p))) {
    return;
  }

  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  const auto populate_upper_triangle =
      _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  if (populate_upper_triangle && p.coords.bin1 <= p.coords.bin2) {
    append(p);
    return;
  }

  if (populate_lower_triangle && p.coords.bin1 >= p.coords.bin2) {
    append(p);
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append_asymmetric(ThinPixel<N> p) {
  assert(!_mirror_pixels);

  if (HICTK_UNLIKELY(!overlaps(p))) {
    return;
  }

  const auto populate_lower_triangle =
      _span == QuerySpan::lower_triangle || _span == QuerySpan::full;
  const auto populate_upper_triangle =
      _span == QuerySpan::upper_triangle || _span == QuerySpan::full;

  if (populate_upper_triangle && p.bin1_id <= p.bin2_id) {
    append(p);
    return;
  }

  if (populate_lower_triangle && p.bin1_id >= p.bin2_id) {
    append(p);
  }
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append(const Pixel<N>& p) {
  assert(_buffer);
  if (_buffer->size() >= _buffer->capacity()) {
    commit_pixels();
  }

  const auto chrom1_id =
      static_cast<std::int32_t>(p.coords.bin1.chrom().id()) - _builder->chrom_id_offset;
  const auto chrom2_id =
      static_cast<std::int32_t>(p.coords.bin2.chrom().id()) - _builder->chrom_id_offset;

  assert(chrom1_id >= 0);
  assert(chrom2_id >= 0);

  _buffer->bin1_id.push_back(static_cast<std::int64_t>(p.coords.bin1.id()));
  _buffer->bin2_id.push_back(static_cast<std::int64_t>(p.coords.bin2.id()));

  _buffer->chrom1_id.push_back(chrom1_id);
  _buffer->start1.push_back(static_cast<std::int32_t>(p.coords.bin1.start()));
  _buffer->end1.push_back(static_cast<std::int32_t>(p.coords.bin1.end()));

  _buffer->chrom2_id.push_back(chrom2_id);
  _buffer->start2.push_back(static_cast<std::int32_t>(p.coords.bin2.start()));
  _buffer->end2.push_back(static_cast<std::int32_t>(p.coords.bin2.end()));

  _buffer->count.push_back(p.count);
}

template <typename PixelIt>
inline void ToDataFrame<PixelIt>::append(const ThinPixel<N>& p) {
  assert(_buffer);
  if (_buffer->size() >= _buffer->capacity()) {
    commit_thin_pixels();
  }

  _buffer->bin1_id.push_back(static_cast<std::int64_t>(p.bin1_id));
  _buffer->bin2_id.push_back(static_cast<std::int64_t>(p.bin2_id));
  _buffer->count.push_back(p.count);
}

template <typename PixelIt>
inline bool ToDataFrame<PixelIt>::overlaps(const ThinPixel<N>& p) const noexcept {
  return overlaps(p.bin1_id, p.bin2_id);
}

template <typename PixelIt>
inline bool ToDataFrame<PixelIt>::overlaps(const Pixel<N>& p) const noexcept {
  return overlaps(p.coords.bin1.id(), p.coords.bin2.id());
}

template <typename PixelIt>
inline bool ToDataFrame<PixelIt>::overlaps(std::uint64_t bin1_id,
                                           std::uint64_t bin2_id) const noexcept {
  if (!_coord1.has_value()) {
    assert(!_coord2.has_value());
    return true;
  }

  // NOLINTBEGIN(*-unchecked-optional-access)
  return bin1_id >= _coord1->bin1.id() && bin1_id <= _coord1->bin2.id() &&
         bin2_id >= _coord2->bin1.id() && bin2_id <= _coord2->bin2.id();
  // NOLINTEND(*-unchecked-optional-access)
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
  assert(_buffer);
  assert(_chunks);

  if (!_buffer->empty()) {
    commit_thin_pixels();
  }

  if (_chunks->empty()) {
    auto result = arrow::Table::MakeEmpty(coo_schema());
    if (!result.ok()) {
      throw std::runtime_error(result.status().ToString());
    }
    return result.MoveValueUnsafe();
  }

  assert(_chunks->bin1_id.size() == _chunks->size());
  assert(_chunks->bin2_id.size() == _chunks->size());
  assert(_chunks->count.size() == _chunks->size());
  auto table =
      arrow::Table::Make(coo_schema(), {std::make_shared<arrow::ChunkedArray>(_chunks->bin1_id),
                                        std::make_shared<arrow::ChunkedArray>(_chunks->bin2_id),
                                        std::make_shared<arrow::ChunkedArray>(_chunks->count)});
  assert(table->ValidateFull().ok());

  _chunks->clear();

  if (_span != QuerySpan::upper_triangle) {
    table = sort_table(table);
  }

  return table;
}

template <typename PixelIt>
inline std::shared_ptr<arrow::Table> ToDataFrame<PixelIt>::make_bg2_table() {
  assert(_buffer);
  assert(_chunks);

  if (!_buffer->empty()) {
    commit_pixels();
  }

  if (_chunks->empty()) {
    auto result = arrow::Table::MakeEmpty(bg2_schema(!_drop_bin_ids));
    if (!result.ok()) {
      throw std::runtime_error(result.status().ToString());
    }
    return result.MoveValueUnsafe();
  }

  auto schema = bg2_schema(!_chunks->bin1_id.empty());
  std::vector<std::shared_ptr<arrow::ChunkedArray>> columns{};
  columns.reserve(static_cast<std::size_t>(schema->num_fields()));

  assert(_chunks->bin1_id.empty() == _chunks->bin2_id.empty());
  if (!_chunks->bin1_id.empty()) {
    columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->bin1_id));
    columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->bin2_id));
  }

  columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->chrom1));
  columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->start1));
  columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->end1));
  columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->chrom2));
  columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->start2));
  columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->end2));
  columns.emplace_back(std::make_shared<arrow::ChunkedArray>(_chunks->count));

  auto table = arrow::Table::Make(std::move(schema), std::move(columns));
  assert(table->ValidateFull().ok());

  _chunks->clear();

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
void ToDataFrame<PixelIt>::commit_thin_pixels() {
  assert(_buffer);
  assert(_builder);
  assert(_chunks);

  if (_buffer->empty()) {
    return;
  }

  append(_builder->int64, _buffer->bin1_id);
  _chunks->bin1_id.emplace_back(finish(_builder->int64));

  append(_builder->int64, _buffer->bin2_id);
  _chunks->bin2_id.emplace_back(finish(_builder->int64));

  append(_builder->count, _buffer->count);
  _chunks->count.emplace_back(finish(_builder->count));

  _buffer->clear();
}

template <typename PixelIt>
void ToDataFrame<PixelIt>::commit_pixels() {
  assert(_buffer);
  assert(_builder);
  assert(_chunks);

  if (_buffer->empty()) {
    return;
  }

  if (_span != QuerySpan::upper_triangle) {
    assert(!_buffer->bin1_id.empty());
    assert(!_buffer->bin2_id.empty());
    append(_builder->int64, _buffer->bin1_id);
    _chunks->bin1_id.emplace_back(finish(_builder->int64));

    append(_builder->int64, _buffer->bin2_id);
    _chunks->bin2_id.emplace_back(finish(_builder->int64));
  }

  append(_builder->chrom, _buffer->chrom1_id);
  _chunks->chrom1.emplace_back(finish(_builder->chrom));

  append(_builder->int32, _buffer->start1);
  _chunks->start1.emplace_back(finish(_builder->int32));

  append(_builder->int32, _buffer->end1);
  _chunks->end1.emplace_back(finish(_builder->int32));

  append(_builder->chrom, _buffer->chrom2_id);
  _chunks->chrom2.emplace_back(finish(_builder->chrom));

  append(_builder->int32, _buffer->start2);
  _chunks->start2.emplace_back(finish(_builder->int32));

  append(_builder->int32, _buffer->end2);
  _chunks->end2.emplace_back(finish(_builder->int32));

  append(_builder->count, _buffer->count);
  _chunks->count.emplace_back(finish(_builder->count));

  _buffer->clear();
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

template <typename PixelIt>
inline QuerySpan ToDataFrame<PixelIt>::fix_query_span(const std::optional<PixelCoordinates>& coord1,
                                                      const std::optional<PixelCoordinates>& coord2,
                                                      QuerySpan requested_span) {
  if (!coord1.has_value()) {
    assert(!coord2.has_value());
    return requested_span;
  }

  // NOLINTNEXTLINE(*-unchecked-optional-access)
  if (coord1->bin2.id() <= coord2->bin1.id()) {
    return QuerySpan::upper_triangle;
  }

  return QuerySpan::full;
}

template <typename PixelIt>
inline std::optional<PixelCoordinates> ToDataFrame<PixelIt>::fix_coordinates(
    std::optional<PixelCoordinates> coord) {
  if (coord.value_or(PixelCoordinates{}) == PixelCoordinates{}) {
    return std::nullopt;
  }

  return coord;
}

}  // namespace hictk::transformers
