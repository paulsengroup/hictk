// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/suppress_warnings.hpp"

namespace hictk::cooler {
template <typename T>
inline Dataset::iterator<T>::iterator(Dataset dset, std::size_t chunk_size, std::size_t h5_offset,
                                      bool init)
    : iterator(std::make_shared<const Dataset>(std::move(dset)), chunk_size, h5_offset, init) {}

template <typename T>
inline Dataset::iterator<T>::iterator(std::shared_ptr<const Dataset> dset, std::size_t chunk_size,
                                      std::size_t h5_offset, bool init)
    // clang-format off
    : _dset(std::move(dset)),
      _h5_chunk_start(h5_offset),
      _h5_offset(h5_offset),
      _chunk_size(chunk_size)
#ifndef NDEBUG
     ,_h5_size(_dset->size())
#endif
// clang-format on
{
  if (_chunk_size == 0) {
    _chunk_size = std::max(std::size_t(2048), _dset->get_chunk_size() / 3);
  }
  if (init) {
    read_chunk_at_offset(_h5_chunk_start);
  }
}

template <typename T>
constexpr bool Dataset::iterator<T>::operator==(const iterator &other) const noexcept {
  return _h5_offset == other._h5_offset;
}
template <typename T>
constexpr bool Dataset::iterator<T>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename T>
constexpr bool Dataset::iterator<T>::operator<(const iterator &other) const noexcept {
  return _h5_offset < other._h5_offset;
}
template <typename T>
constexpr bool Dataset::iterator<T>::operator<=(const iterator &other) const noexcept {
  return _h5_offset <= other._h5_offset;
}

template <typename T>
constexpr bool Dataset::iterator<T>::operator>(const iterator &other) const noexcept {
  return _h5_offset > other._h5_offset;
}
template <typename T>
constexpr bool Dataset::iterator<T>::operator>=(const iterator &other) const noexcept {
  return _h5_offset >= other._h5_offset;
}

template <typename T>
inline auto Dataset::iterator<T>::operator*() const -> value_type {
  switch (underlying_buff_status()) {
    case OVERLAPPING:
      break;
    case UNINITIALIZED:
      [[fallthrough]];
    case DOWNSTEAM:
      // Read first chunk
      read_chunk_at_offset(_h5_offset);
      break;
    case UPSTREAM:
      // Iterator was decremented one or more times since the last dereference, thus we assume the
      // iterator is being used to traverse the dataset backward
      _h5_chunk_start = _h5_offset - (std::min)(_buff->size() - 1, _h5_offset);
      read_chunk_at_offset(_h5_chunk_start);
  }

  assert(_buff);
  assert(_h5_offset < _h5_size);
  assert(_h5_chunk_start <= _h5_offset);
  assert(_h5_offset - _h5_chunk_start < _buff->size());
  return (*_buff)[_h5_offset - _h5_chunk_start];
}

template <typename T>
inline auto Dataset::iterator<T>::operator[](std::size_t i) const -> value_type {
  return *(*this + i);
}

template <typename T>
inline auto Dataset::iterator<T>::operator++() -> iterator & {
  return (*this) += 1;
}

template <typename T>
inline auto Dataset::iterator<T>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  if (_h5_offset > _h5_chunk_start + _chunk_size) {
    read_chunk_at_offset(_h5_offset);
  }
  return it;
}

template <typename T>
inline auto Dataset::iterator<T>::operator+=(difference_type i) -> iterator & {
  if (HICTK_UNLIKELY(i < 0)) {
    return *this -= -i;
  }

  HICTK_DISABLE_WARNING_PUSH
  HICTK_DISABLE_WARNING_SIGN_COMPARE
  HICTK_DISABLE_WARNING_SIGN_CONVERSION
  HICTK_DISABLE_WARNING_CONVERSION
  assert(_h5_offset + i <= _h5_size);
  _h5_offset += i;
  HICTK_DISABLE_WARNING_POP
  return *this;
}

template <typename T>
inline auto Dataset::iterator<T>::operator+(difference_type i) const -> iterator {
  if (HICTK_UNLIKELY(i < 0)) {
    return *this - -i;
  }

  assert(_buff);
  HICTK_DISABLE_WARNING_PUSH
  HICTK_DISABLE_WARNING_SIGN_COMPARE
  HICTK_DISABLE_WARNING_SIGN_CONVERSION
  HICTK_DISABLE_WARNING_CONVERSION
  const auto new_offset = _h5_offset + i;
  assert(new_offset <= _h5_size);

  if (!_buff || _h5_chunk_start + _buff->size() < new_offset) {
    return iterator(*_dset, _chunk_size, new_offset);
  }
  HICTK_DISABLE_WARNING_POP

  auto it = *this;
  return it += i;
}

template <typename T>
inline auto Dataset::iterator<T>::operator--() -> iterator & {
  assert(_h5_offset != 0);
  return (*this) -= 1;
}

template <typename T>
inline auto Dataset::iterator<T>::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  if (_h5_offset < _h5_chunk_start) {
    read_chunk_at_offset(_h5_offset - (std::min)(_chunk_size - 1, _h5_offset));
  }
  return it;
}

template <typename T>
inline auto Dataset::iterator<T>::operator-=(difference_type i) -> iterator & {
  if (HICTK_UNLIKELY(i < 0)) {
    return *this += -i;
  }

  HICTK_DISABLE_WARNING_PUSH
  HICTK_DISABLE_WARNING_SIGN_COMPARE
  HICTK_DISABLE_WARNING_SIGN_CONVERSION
  HICTK_DISABLE_WARNING_CONVERSION
  assert(_h5_offset >= i);
  _h5_offset -= i;
  HICTK_DISABLE_WARNING_POP
  return *this;
}

template <typename T>
inline auto Dataset::iterator<T>::operator-(difference_type i) const -> iterator {
  if (HICTK_UNLIKELY(i < 0)) {
    return *this + -i;
  }

  HICTK_DISABLE_WARNING_PUSH
  HICTK_DISABLE_WARNING_SIGN_COMPARE
  HICTK_DISABLE_WARNING_SIGN_CONVERSION
  HICTK_DISABLE_WARNING_CONVERSION
  assert(_h5_offset >= i);
  const auto new_offset = _h5_offset - i;
  if (new_offset >= _h5_chunk_start) {
    auto it = *this;
    return it -= i;
  }
  HICTK_DISABLE_WARNING_POP

  return iterator(*_dset, _chunk_size, new_offset);
}

template <typename T>
inline auto Dataset::iterator<T>::operator-(const iterator &other) const -> difference_type {
  return static_cast<difference_type>(_h5_offset) - static_cast<difference_type>(other._h5_offset);
}

template <typename T>
inline auto Dataset::iterator<T>::seek(std::size_t offset) -> iterator<T> & {
  assert(offset < _h5_size);
  if (offset >= h5_offset()) {
    return *this += static_cast<std::ptrdiff_t>(offset - h5_offset());
  }
  return *this -= static_cast<std::ptrdiff_t>(h5_offset() - offset);
}

template <typename T>
constexpr std::uint64_t Dataset::iterator<T>::h5_offset() const noexcept {
  return _h5_offset;
}

template <typename T>
constexpr std::size_t Dataset::iterator<T>::underlying_buff_capacity() const noexcept {
  return _chunk_size;
}

template <typename T>
constexpr std::size_t Dataset::iterator<T>::lower_bound() const noexcept {
  return _h5_chunk_start;
}

template <typename T>
constexpr std::size_t Dataset::iterator<T>::upper_bound() const noexcept {
  if (_buff) {
    return _h5_chunk_start + _buff->size();
  }
  return _h5_chunk_start + _chunk_size;
}

template <typename T>
constexpr auto Dataset::iterator<T>::underlying_buff_status() const noexcept -> OverlapStatus {
  if (!_buff) {
    return UNINITIALIZED;
  }

  if (_h5_offset >= upper_bound()) {
    return DOWNSTEAM;
  }

  if (_h5_offset - lower_bound() >= _buff->size()) {
    return UPSTREAM;
  }

  return OVERLAPPING;
}

template <typename T>
constexpr std::size_t Dataset::iterator<T>::underlying_buff_num_available_rev() const noexcept {
  if (underlying_buff_status() != OVERLAPPING) {
    return 0;
  }
  return _h5_offset - lower_bound();
}

template <typename T>
constexpr std::size_t Dataset::iterator<T>::underlying_buff_num_available_fwd() const noexcept {
  if (underlying_buff_status() != OVERLAPPING) {
    return 0;
  }
  return upper_bound() - _h5_offset;
}

template <typename T>
constexpr const Dataset &Dataset::iterator<T>::dataset() const noexcept {
  assert(_dset);
  return *_dset;
}

template <typename T>
inline void Dataset::iterator<T>::read_chunk_at_offset(std::size_t new_offset) const {
  assert(_dset);

  const auto dset_size = dataset().size();

  if (new_offset == dset_size) {
    _buff = nullptr;
    _h5_chunk_start = dset_size;
    return;
  }

  if (!_buff || _buff.use_count() > 1) {
    //  This should be fine, as copying Dataset::iterator is not thread-safe anyway
    _buff = std::make_shared<std::vector<T>>(_chunk_size);
  }

  const auto buff_size = (std::min)(_chunk_size, dataset().size() - new_offset);
  _buff->resize(buff_size);
  _dset->read(*_buff, buff_size, new_offset);

  _h5_chunk_start = new_offset;
}

template <typename T>
inline auto Dataset::iterator<T>::make_end_iterator(Dataset dset, std::size_t chunk_size)
    -> iterator {
  return iterator::make_end_iterator(std::make_shared<const Dataset>(std::move(dset)), chunk_size);
}

template <typename T>
inline auto Dataset::iterator<T>::make_end_iterator(std::shared_ptr<const Dataset> dset,
                                                    std::size_t chunk_size) -> iterator {
  iterator it{};
  it._buff = nullptr;
  it._dset = std::move(dset);
  it._h5_offset = it._dset->size();
  it._chunk_size =
      chunk_size == 0 ? std::max(std::size_t(2048), it._dset->get_chunk_size() / 3) : chunk_size;
#ifndef NDEBUG
  it._h5_size = it._h5_offset;
#endif
  it._h5_chunk_start = it._h5_offset;

  return it;
}

}  // namespace hictk::cooler
