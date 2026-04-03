// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/pixel_selector.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/pixel.hpp"

namespace hictk::cooler {

PixelSelector::PixelSelector(const Index &index, const Dataset &pixels_bin1_id,
                             const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                             std::shared_ptr<const balancing::Weights> weights,
                             bool symmetric_upper_) noexcept
    : PixelSelector(index.bins_ptr(), pixels_bin1_id, pixels_bin2_id, pixels_count,
                    std::move(weights), symmetric_upper_) {}

PixelSelector::PixelSelector(std::shared_ptr<const BinTable> bins, const Dataset &pixels_bin1_id,
                             const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                             std::shared_ptr<const balancing::Weights> weights,
                             bool symmetric_upper_) noexcept
    : _bins(std::move(bins)),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count),
      _weights(std::move(weights)),
      _symmetric_upper(symmetric_upper_) {
  assert(_bins);
  assert(_weights);
}

PixelSelector::PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                             const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                             const PixelCoordinates &coords,
                             std::shared_ptr<const balancing::Weights> weights,
                             bool symmetric_upper_)
    : PixelSelector(std::move(index), pixels_bin1_id, pixels_bin2_id, pixels_count, coords, coords,
                    std::move(weights), symmetric_upper_) {}

PixelSelector::PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                             const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                             PixelCoordinates coord1, PixelCoordinates coord2,
                             std::shared_ptr<const balancing::Weights> weights,
                             bool symmetric_upper_)
    : _coord1(std::move(coord1)),
      _coord2(std::move(coord2)),
      _index(std::move(index)),
      _bins(_index->bins_ptr()),
      _pixels_bin1_id(&pixels_bin1_id),
      _pixels_bin2_id(&pixels_bin2_id),
      _pixels_count(&pixels_count),
      _weights(std::move(weights)),
      _symmetric_upper(symmetric_upper_) {
  assert(_index);
  assert(_weights);

  if (!_symmetric_upper) {
    assert(_coord1.empty() && _coord2.empty());
  }

  const auto query_is_cis = _coord1.bin1.chrom() == _coord2.bin1.chrom();
  if ((!query_is_cis && _coord1.bin1 > _coord2.bin1) ||
      (query_is_cis && _coord1.bin1.start() > _coord2.bin1.start())) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("query {}:{}-{}; {}:{}-{}; overlaps with the lower-triangle of the matrix"),
        _coord1.bin1.chrom().name(), _coord1.bin1.start(), _coord1.bin2.end(),
        _coord2.bin1.chrom().name(), _coord2.bin1.start(), _coord2.bin2.end()));
  }
}

bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  // clang-format off
  return begin<int>() == other.begin<int>() &&
         end<int>() == other.end<int>() &&
         _weights == other._weights &&
         _symmetric_upper == other._symmetric_upper;
  // clang-format on
}

bool PixelSelector::operator!=(const PixelSelector &other) const noexcept {
  return !(*this == other);
}

bool PixelSelector::empty() const { return begin<double>() == end<double>(); }

const PixelCoordinates &PixelSelector::coord1() const noexcept { return _coord1; }

const PixelCoordinates &PixelSelector::coord2() const noexcept { return _coord2; }

std::uint64_t PixelSelector::size(bool upper_triangle) const {
  if (!_coord1) {
    assert(!_coord2);
    const auto n = bins().size();
    if (upper_triangle && _symmetric_upper) {
      return (n * (n + 1)) / 2;
    }
    return n * n;
  }

  if (bins().type() != BinTable::Type::fixed) {
    throw std::runtime_error(
        "computing the number of pixels overlapping a query over a matrix with variable bin size "
        "is not supported.");
  }

  assert(_symmetric_upper);

  if (_coord1.bin1.chrom() != _coord2.bin1.chrom()) {
    const auto height = _coord1.bin2.id() - _coord1.bin1.id() + 1;
    const auto width = _coord2.bin2.id() - _coord2.bin1.id() + 1;
    return height * width;
  }

  const auto start1 = _coord1.bin1.start();
  const auto start2 = _coord2.bin1.start();

  const auto end1 = ((_coord1.bin2.rel_id() + 1) * bins().resolution()) + 1;
  const auto end2 = ((_coord2.bin2.rel_id() + 1) * bins().resolution()) + 1;

  return area(start1, end1, start2, end2, bins().resolution(), upper_triangle);
}

const BinTable &PixelSelector::bins() const noexcept { return *bins_ptr(); }

std::shared_ptr<const BinTable> PixelSelector::bins_ptr() const noexcept { return _bins; }

PixelSelector PixelSelector::fetch(PixelCoordinates coord1, PixelCoordinates coord2) const {
  return {_index,         *_pixels_bin1_id,  *_pixels_bin2_id,
          *_pixels_count, std::move(coord1), std::move(coord2),
          _weights,       _symmetric_upper};
}

const balancing::Weights &PixelSelector::weights() const noexcept {
  assert(_weights);
  return *_weights;
}

bool PixelSelector::is_symmetric_upper() const noexcept { return _symmetric_upper; }

}  // namespace hictk::cooler
