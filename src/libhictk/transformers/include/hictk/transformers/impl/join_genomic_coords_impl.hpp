// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename PixelIt>
inline JoinGenomicCoords<PixelIt>::JoinGenomicCoords(PixelIt first, PixelIt last,
                                                     std::shared_ptr<const BinTable> bins)
    : _first(std::move(first)), _last(std::move(last)), _bins(std::move(bins)) {}

template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::begin() const -> iterator {
  return iterator{_first, _bins};
}
template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::end() const -> iterator {
  return iterator::at_end(_last, _bins);
}

template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::cbegin() const -> iterator {
  return begin();
}
template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::cend() const -> iterator {
  return end();
}

template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::read_all() const -> std::vector<Pixel<N>> {
  assert(_bins);
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel<N>> buff{};
  std::transform(_first, _last, std::back_inserter(buff), [&](const ThinPixel<N>& p) {
    return Pixel<N>{{_bins->at(p.bin1_id), _bins->at(p.bin2_id)}, p.count};
  });
  return buff;
}

template <typename PixelIt>
inline JoinGenomicCoords<PixelIt>::iterator::iterator(PixelIt it,
                                                      std::shared_ptr<const BinTable> bins)
    : _it(std::move(it)), _bins(std::move(bins)) {}

template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::iterator::at_end(PixelIt it,
                                                         std::shared_ptr<const BinTable> bins)
    -> iterator {
  iterator it_{};
  it_._it = std::move(it);
  it_._bins = std::move(bins);
  return it_;
}

template <typename PixelIt>
inline bool JoinGenomicCoords<PixelIt>::iterator::operator==(const iterator& other) const noexcept {
  return _it == other._it;
}
template <typename PixelIt>
inline bool JoinGenomicCoords<PixelIt>::iterator::operator!=(const iterator& other) const noexcept {
  return !(*this == other);
}

template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::iterator::operator*() const -> const_reference {
  assert(_bins);
  _value = Pixel<N>{_bins->at(_it->bin1_id), _bins->at(_it->bin2_id), _it->count};
  return _value;
}
template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::iterator::operator->() const -> const_pointer {
  return &(**this);
}

template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::iterator::operator++() -> iterator& {
  ++_it;
  return *this;
}
template <typename PixelIt>
inline auto JoinGenomicCoords<PixelIt>::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++_it;
  return it;
}
}  // namespace hictk::transformers
