// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/common.hpp"

namespace hictk::cooler {

namespace internal {
template <typename T>
inline COWChunk<T>::COWChunk(std::size_t start_, SharedBufferT data_,
                             std::size_t capacity_) noexcept
    : _buff(std::move(data_)), _start(start_) {
  reserve(capacity_);
}

template <typename T>
inline COWChunk<T>::COWChunk(std::size_t start_, BufferT data_, std::size_t capacity_)
    : COWChunk(start_, data_.empty() ? nullptr : std::make_shared<BufferT>(std::move(data_)),
               capacity_) {}

template <typename T>
constexpr std::size_t COWChunk<T>::COWChunk::id() const noexcept {
  if (capacity() == 0) {
    return 0;
  }

  return _start / capacity();
}

template <typename T>
constexpr std::size_t COWChunk<T>::COWChunk::start() const noexcept {
  return _start;
}

template <typename T>
inline std::size_t COWChunk<T>::COWChunk::end() const noexcept {
  return _start + size();
}

template <typename T>
inline std::size_t COWChunk<T>::COWChunk::capacity() const noexcept {
  return !_buff ? 0 : _buff->capacity();
}

template <typename T>
inline std::size_t COWChunk<T>::COWChunk::size() const noexcept {
  return !_buff ? 0 : _buff->size();
}

template <typename T>
inline bool COWChunk<T>::COWChunk::empty() const noexcept {
  return size() == 0;
}

template <typename T>
inline std::size_t COWChunk<T>::COWChunk::use_count() const noexcept {
  return static_cast<std::size_t>(std::max(0L, _buff.use_count()));
}

template <typename T>
inline auto COWChunk<T>::operator()() const noexcept -> const BufferT & {
  if (!!_buff) {
    return *_buff;
  }

  return _empty_buffer;
}

template <typename T>
inline auto COWChunk<T>::operator()() noexcept -> BufferT & {
  assert(!!_buff);
  return *_buff;  // NOLINT
}

template <typename T>
inline auto COWChunk<T>::operator()(std::size_t i) const noexcept -> std::optional<T> {
  if (i >= start() && i < end()) {
    assert(!!_buff);
    return (*_buff)[i - start()];  // NOLINT
  }

  return {};
}

template <typename T>
inline auto COWChunk<T>::operator[](std::size_t i) const noexcept -> T {
  assert(!!_buff);
  assert(i >= start());
  assert(i < end());
  return (*_buff)[i - start()];  // NOLINT
}

template <typename T>
inline void COWChunk<T>::update(std::size_t start_) noexcept {
  _start = start_;
}

template <typename T>
inline void COWChunk<T>::update(std::size_t start_, SharedBufferT data_) {
  if (data_ && data_->capacity() != capacity()) {
    throw std::runtime_error(
        "COWChunk<T>::update(): incoming data has a different size then the current buffer! "
        "Hint: call resize() or reserve() before calling update()");
  }
  update(start_);
  _buff = std::move(data_);
}

template <typename T>
inline void COWChunk<T>::update(std::size_t start_, BufferT data_) {
  update(start_);
  if (data_.empty()) {
    reset_buffer();
    return;
  }

  if (data_.size() > capacity()) {
    throw std::runtime_error(
        "COWChunk<T>::update(): incoming data is larger than the available space! "
        "Hint: call resize() or reserve() before calling update()");
  }

  if (!_buff || use_count() > 1) {
    data_.reserve(capacity());
    update(start_, std::make_shared<BufferT>(std::move(data_)));
    return;
  }

  _buff->resize(data_.size());
  std::move(data_.begin(), data_.end(), _buff->begin());
}

template <typename T>
inline void COWChunk<T>::resize(std::size_t new_size, bool shrink_to_fit) {
  if (new_size == size() && !shrink_to_fit) {
    return;
  }
  if (new_size == 0) {
    _buff = nullptr;
  }

  if (!_buff) {
    // need to allocate a new buffer
    _buff = std::make_shared<BufferT>(new_size);
  } else if (use_count() > 1) {
    // need to allocate a new buffer and copy the old data
    auto new_buff = std::make_shared<BufferT>(new_size);
    const auto end_offset = static_cast<std::ptrdiff_t>(std::min(_buff->size(), new_buff->size()));
    std::move(_buff->begin(), _buff->begin() + end_offset, new_buff->begin());
    _buff = std::move(new_buff);
  } else {
    // we can resize the current buffer
    _buff->resize(new_size);
    if (shrink_to_fit) {
      _buff->shrink_to_fit();
    }
  }
}

template <typename T>
inline void COWChunk<T>::reserve(std::size_t new_capacity) {
  if (!_buff) {
    _buff = std::make_shared<BufferT>();
  }
  _buff->reserve(new_capacity);
}

template <typename T>
inline void COWChunk<T>::reset_buffer() noexcept {
  _buff = nullptr;
}

}  // namespace internal

template <typename T>
inline Dataset::iterator<T>::iterator(Dataset dset, std::optional<std::ptrdiff_t> chunk_size_,
                                      std::size_t h5_offset, bool init)
    : iterator(std::make_shared<const Dataset>(std::move(dset)), chunk_size_, h5_offset, init) {}

template <typename T>
inline Dataset::iterator<T>::iterator(std::shared_ptr<const Dataset> dset,
                                      std::optional<std::ptrdiff_t> chunk_size_,
                                      std::size_t h5_offset, bool init)

    : _buffer(h5_offset, nullptr),
      _dset(std::move(dset)),
      _chunk_size(compute_chunk_size(_dset, chunk_size_)),
      _h5_offset(h5_offset),
      _h5_size(!!_dset ? _dset->size() : 0) {
  if (init) {
    const auto read_forward = chunk_size_.value_or(0) >= 0;
    read_chunk_at_offset(_h5_offset, read_forward);
    assert(_h5_offset >= _buffer.start());
    assert(_buffer.empty() || _h5_offset <= _buffer.end());
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
  bound_check();
  if (HICTK_UNLIKELY(buffer_is_outdated())) {
    read_chunk_at_offset(_h5_offset, _h5_offset >= _buffer.end());
  }
  assert(_buffer.start() <= _h5_offset);
  assert(_h5_offset < _buffer.end());

  return _buffer[_h5_offset];
}

template <typename T>
inline auto Dataset::iterator<T>::operator[](std::size_t i) const -> value_type {
  return *(*this + i);
}

template <typename T>
inline auto Dataset::iterator<T>::operator++() -> iterator & {
  *this += 1;
  if (_h5_offset >= _buffer.end()) {
    read_chunk_at_offset(_h5_offset, true);
  }
  return *this;
}

template <typename T>
inline auto Dataset::iterator<T>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename T>
inline auto Dataset::iterator<T>::operator+=(difference_type i) -> iterator & {
  if (HICTK_UNLIKELY(i < 0)) {
    return *this -= -i;
  }

  bound_check(i, true);
  _h5_offset += static_cast<std::size_t>(i);
  return *this;
}

template <typename T>
inline auto Dataset::iterator<T>::operator+(difference_type i) const -> iterator {
  if (HICTK_UNLIKELY(i < 0)) {
    return *this - -i;
  }

  bound_check(i, true);
  const auto new_offset = _h5_offset + static_cast<std::size_t>(i);

  if (!_buffer.empty() && _buffer.end() < new_offset) {
    return iterator{*_dset, static_cast<std::ptrdiff_t>(_chunk_size), new_offset};
  }

  auto it = *this;
  return it += i;
}

template <typename T>
inline auto Dataset::iterator<T>::operator--() -> iterator & {
  *this -= 1;
  if (_h5_offset < _buffer.start()) {
    read_chunk_at_offset(_h5_offset, false);
  }
  return *this;
}

template <typename T>
inline auto Dataset::iterator<T>::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  return it;
}

template <typename T>
inline auto Dataset::iterator<T>::operator-=(difference_type i) -> iterator & {
  if (HICTK_UNLIKELY(i < 0)) {
    return *this += -i;
  }

  bound_check(-i);
  _h5_offset -= static_cast<std::size_t>(i);
  return *this;
}

template <typename T>
inline auto Dataset::iterator<T>::operator-(difference_type i) const -> iterator {
  if (HICTK_UNLIKELY(i < 0)) {
    return *this + -i;
  }

  bound_check(-i);
  const auto new_offset = _h5_offset - static_cast<std::size_t>(i);
  if (new_offset >= _buffer.start()) {
    auto it = *this;
    return it -= i;
  }

  return iterator{*_dset, static_cast<std::ptrdiff_t>(_chunk_size), new_offset};
}

template <typename T>
inline auto Dataset::iterator<T>::operator-(const iterator &other) const -> difference_type {
  return static_cast<difference_type>(_h5_offset) - static_cast<difference_type>(other._h5_offset);
}

template <typename T>
template <typename I>
inline auto Dataset::iterator<T>::seek(I offset) -> iterator & {
  static_assert(std::is_integral_v<I>);
  const auto rel_offset =
      conditional_static_cast<std::ptrdiff_t>(offset) - static_cast<std::ptrdiff_t>(_h5_offset);
  bound_check(rel_offset, true);

  return *this += rel_offset;
}

template <typename T>
constexpr std::size_t Dataset::iterator<T>::h5_offset() const noexcept {
  return _h5_offset;
}

template <typename T>
inline auto Dataset::iterator<T>::buffer() const -> const internal::COWChunk<T> & {
  if (HICTK_UNLIKELY(_buffer.empty() && _chunk_size != 0)) {
    read_chunk_at_offset(_buffer.start(), true);
  } else if (HICTK_UNLIKELY(_h5_offset < _buffer.start())) {
    read_chunk_at_offset(_h5_offset, false);
  } else if (HICTK_UNLIKELY(_h5_offset >= _buffer.end()) && _h5_offset != _h5_size) {
    read_chunk_at_offset(_h5_offset, true);
  }
  return _buffer;
}

template <typename T>
constexpr std::size_t Dataset::iterator<T>::chunk_size() const noexcept {
  return static_cast<std::size_t>(_chunk_size);
}

template <typename T>
constexpr const Dataset &Dataset::iterator<T>::dataset() const noexcept {
  assert(!!_dset);
  return *_dset;
}

template <typename T>
inline void Dataset::iterator<T>::read_chunk_at_offset(std::size_t new_offset, bool forward) const {
  assert(!!_dset);

  if (_chunk_size == 0) {
    _buffer.update(new_offset, nullptr);
    return;
  }

  if (new_offset >= _h5_size && forward) {
    // we are at/past the end
    _buffer.update(_h5_size, nullptr);
    return;
  }

  assert(_chunk_size != 0);
  const auto start_offset = (new_offset / _chunk_size) * _chunk_size;

  std::size_t size = 0;
  if (start_offset < _h5_size) {
    size = std::min(static_cast<std::size_t>(_chunk_size), _h5_size - start_offset);
  }

  assert(new_offset >= start_offset);
  if (new_offset < _h5_size) {
    assert(new_offset < start_offset + size);
  } else {
    assert(start_offset + size == _h5_size);
  }

  if (_buffer.start() == start_offset && _buffer.size() >= size) {
    return;
  }

  if (_buffer.use_count() > 1) {
    _buffer.reset_buffer();
  }
  _buffer.resize(size);
  _dset->read(_buffer(), _buffer.size(), start_offset);
  _buffer.update(start_offset);
}

template <typename T>
inline bool Dataset::iterator<T>::buffer_is_outdated() const noexcept {
  if (_h5_offset >= _h5_size) {
    return false;
  }

  if (_buffer.empty()) {
    return true;
  }

  if (_h5_offset < _buffer.start()) {
    return true;
  }

  if (_h5_offset >= _buffer.end()) {
    return true;
  }

  return false;
}

template <typename T>
inline auto Dataset::iterator<T>::make_end_iterator(Dataset dset,
                                                    std::optional<std::ptrdiff_t> chunk_size_)
    -> iterator {
  return iterator::make_end_iterator(std::make_shared<const Dataset>(std::move(dset)), chunk_size_);
}

template <typename T>
inline auto Dataset::iterator<T>::make_end_iterator(std::shared_ptr<const Dataset> dset,
                                                    std::optional<std::ptrdiff_t> chunk_size_)
    -> iterator {
  const auto offset = dset->size();
  return iterator{std::move(dset), chunk_size_, offset, chunk_size_ != 0};
}

template <typename T>
inline std::uint32_t Dataset::iterator<T>::compute_chunk_size(
    const std::shared_ptr<const Dataset> &dset, std::optional<std::ptrdiff_t> chunk_size_) {
  if (!dset) {
    return 0;
  }

  const auto size = chunk_size_.value_or(static_cast<std::ptrdiff_t>(dset->get_chunk_size()));
  constexpr auto max_size = static_cast<std::ptrdiff_t>(std::numeric_limits<std::uint32_t>::max());
  return static_cast<std::uint32_t>(std::min(std::abs(size), max_size));
}

template <typename T>
inline void Dataset::iterator<T>::bound_check([[maybe_unused]] std::ptrdiff_t i,
                                              [[maybe_unused]] bool close_interval) const noexcept {
#ifndef NDEBUG
  if (i < 0) {
    assert(static_cast<std::ptrdiff_t>(_h5_offset) >= -i);
    return;
  }

  i += static_cast<std::ptrdiff_t>(_h5_offset);
  if (close_interval) {
    assert(i <= static_cast<std::ptrdiff_t>(_h5_size));
  } else {
    assert(i < static_cast<std::ptrdiff_t>(_h5_size));
  }
#endif
}

}  // namespace hictk::cooler
