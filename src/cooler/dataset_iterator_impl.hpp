// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Selection.hpp>
#include <memory>
#include <string>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/type_pretty_printer.hpp"

namespace hictk {

template <typename T, std::size_t CHUNK_SIZE>
inline Dataset::iterator<T, CHUNK_SIZE>::iterator(const Dataset &dset, std::size_t h5_offset,
                                                  bool init)
    // clang-format off
    : _dset(&dset),
      _h5_chunk_start(h5_offset),
      _h5_offset(h5_offset)
#ifndef NDEBUG
     ,_h5_size(dset.size())
#endif
// clang-format on
{
  if (init) {
    this->read_chunk_at_offset(this->_h5_chunk_start);
  }
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr bool Dataset::iterator<T, CHUNK_SIZE>::operator==(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset == other._h5_offset;
}
template <typename T, std::size_t CHUNK_SIZE>
constexpr bool Dataset::iterator<T, CHUNK_SIZE>::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr bool Dataset::iterator<T, CHUNK_SIZE>::operator<(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset < other._h5_offset;
}
template <typename T, std::size_t CHUNK_SIZE>
constexpr bool Dataset::iterator<T, CHUNK_SIZE>::operator<=(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset <= other._h5_offset;
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr bool Dataset::iterator<T, CHUNK_SIZE>::operator>(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset > other._h5_offset;
}
template <typename T, std::size_t CHUNK_SIZE>
constexpr bool Dataset::iterator<T, CHUNK_SIZE>::operator>=(const iterator &other) const noexcept {
  assert(this->_dset == other._dset);
  return this->_h5_offset >= other._h5_offset;
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator*() const -> value_type {
  switch (this->underlying_buff_status()) {
    case OVERLAPPING:
      break;
    case UNINITIALIZED:
      [[fallthrough]];
    case DOWNSTEAM:
      // Read first chunk
      this->read_chunk_at_offset(this->_h5_offset);
      break;
    case UPSTREAM:
      // Iterator was decremented one or more times since the last dereference, thus we assume the
      // iterator is being used to traverse the dataset backward
      this->_h5_chunk_start =
          this->_h5_offset - (std::min)(this->_buff->size() - 1, this->_h5_offset);
      this->read_chunk_at_offset(this->_h5_chunk_start);
  }

  assert(this->_buff);
  assert(this->_dset);
  assert(this->_h5_offset < this->_h5_size);
  assert(this->_h5_chunk_start <= this->_h5_offset);
  assert(this->_h5_offset - this->_h5_chunk_start < this->_buff->size());
  return (*this->_buff)[this->_h5_offset - this->_h5_chunk_start];
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator[](std::size_t i) const -> value_type {
  return *(*this + i);
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator++() -> iterator & {
  return (*this) += 1;
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  if (this->_h5_offset > this->_h5_chunk_start + CHUNK_SIZE) {
    this->read_chunk_at_offset(this->_h5_offset);
  }
  return it;
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator+=(std::size_t i) -> iterator & {
  assert(this->_dset);
  assert(this->_h5_offset + i <= this->_h5_size);
  this->_h5_offset += i;
  return *this;
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator+(std::size_t i) const -> iterator {
  assert(this->_dset);
  assert(this->_buff);
  const auto new_offset = this->_h5_offset + i;
  assert(new_offset <= this->_h5_size);

  if (!this->_buff || this->_h5_chunk_start + this->_buff->size() < new_offset) {
    return iterator(*this->_dset, new_offset);
  }

  auto it = *this;
  return it += i;
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator--() -> iterator & {
  assert(this->_h5_offset != 0);
  return (*this) -= 1;
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  if (this->_h5_offset < this->_h5_chunk_start) {
    this->read_chunk_at_offset(this->_h5_offset - (std::min)(CHUNK_SIZE - 1, this->_h5_offset));
  }
  return it;
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator-=(std::size_t i) -> iterator & {
  assert(this->_h5_offset >= i);
  this->_h5_offset -= i;
  return *this;
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator-(std::size_t i) const -> iterator {
  assert(this->_h5_offset >= i);
  const auto new_offset = this->_h5_offset - i;
  if (new_offset >= this->_h5_chunk_start) {
    auto it = *this;
    return it -= i;
  }

  assert(this->_dset);
  return iterator(*this->_dset, new_offset);
}

template <typename T, std::size_t CHUNK_SIZE>
inline auto Dataset::iterator<T, CHUNK_SIZE>::operator-(const iterator &other) const
    -> difference_type {
  return static_cast<difference_type>(this->_h5_offset) -
         static_cast<difference_type>(other._h5_offset);
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr std::uint64_t Dataset::iterator<T, CHUNK_SIZE>::h5_offset() const noexcept {
  return this->_h5_offset;
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr std::size_t Dataset::iterator<T, CHUNK_SIZE>::underlying_buff_capacity() const noexcept {
  return CHUNK_SIZE;
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr std::size_t Dataset::iterator<T, CHUNK_SIZE>::lower_bound() const noexcept {
  return this->_h5_chunk_start;
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr std::size_t Dataset::iterator<T, CHUNK_SIZE>::upper_bound() const noexcept {
  if (this->_buff) {
    return this->_h5_chunk_start + this->_buff->size();
  }
  return this->_h5_chunk_start + CHUNK_SIZE;
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr auto Dataset::iterator<T, CHUNK_SIZE>::underlying_buff_status() const noexcept
    -> OverlapStatus {
  if (!this->_buff) {
    return UNINITIALIZED;
  }

  if (this->_h5_offset >= this->upper_bound()) {
    return DOWNSTEAM;
  }

  if (this->_h5_offset - this->lower_bound() >= this->_buff->size()) {
    return UPSTREAM;
  }

  return OVERLAPPING;
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr std::size_t Dataset::iterator<T, CHUNK_SIZE>::underlying_buff_num_available_rev()
    const noexcept {
  if (this->underlying_buff_status() != OVERLAPPING) {
    return 0;
  }
  return this->_h5_offset - this->lower_bound();
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr std::size_t Dataset::iterator<T, CHUNK_SIZE>::underlying_buff_num_available_fwd()
    const noexcept {
  if (this->underlying_buff_status() != OVERLAPPING) {
    return 0;
  }
  return this->upper_bound() - this->_h5_offset;
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr const Dataset &Dataset::iterator<T, CHUNK_SIZE>::dataset() const noexcept {
  return *this->_dset;
}

template <typename T, std::size_t CHUNK_SIZE>
inline void Dataset::iterator<T, CHUNK_SIZE>::read_chunk_at_offset(std::size_t new_offset) const {
  assert(this->_dset);

  const auto dset_size = this->dataset().size();

  if (new_offset == dset_size) {
    this->_buff = nullptr;
    this->_h5_chunk_start = dset_size;
    return;
  }

  if (!this->_buff || this->_buff.use_count() > 1) {
    //  This should be fine, as copying Dataset::iterator is not thread-safe anyway
    this->_buff = std::make_shared<std::vector<T>>(CHUNK_SIZE);
  }

  const auto buff_size = (std::min)(CHUNK_SIZE, this->dataset().size() - new_offset);
  this->_buff->resize(buff_size);
  this->_dset->read(*this->_buff, buff_size, new_offset);

  this->_h5_chunk_start = new_offset;
}

template <typename T, std::size_t CHUNK_SIZE>
constexpr auto Dataset::iterator<T, CHUNK_SIZE>::make_end_iterator(const Dataset &dset)
    -> iterator {
  iterator it{};
  it._buff = nullptr;
  it._dset = &dset;
  it._h5_offset = dset.size();
#ifndef NDEBUG
  it._h5_size = it._h5_offset;
#endif
  it._h5_chunk_start = it._h5_offset;

  return it;
}

}  // namespace hictk
