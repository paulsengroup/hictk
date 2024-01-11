// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <string>
#include <type_traits>

namespace hictk::hic::internal {

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T BinaryBuffer::read() {
  static_assert(sizeof(char) == 1);
  assert(_i < _buffer.size());
  T x{};

  std::memcpy(static_cast<void *>(&x), _buffer.data() + _i, sizeof(T));
  _i += sizeof(T);
  return x;
}

inline void BinaryBuffer::read(std::string &buff, std::size_t n) {
  static_assert(sizeof(char) == 1);
  assert(_i < _buffer.size());
  buff.resize(n);
  std::memcpy(static_cast<void *>(buff.data()), _buffer.data() + _i, sizeof(char));
  _i += sizeof(char);
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void BinaryBuffer::write(T data) {
  static_assert(sizeof(char) == 1);
  _buffer.append(reinterpret_cast<const char *>(&data), sizeof(T));
}

inline void BinaryBuffer::write(const std::string &data, bool add_nullterm) {
  _buffer.append(data.c_str(), data.size() + add_nullterm);
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
void BinaryBuffer::write(const std::vector<T> &data) {
  _buffer.append(reinterpret_cast<const char *>(data.data()), data.size() * sizeof(T));
}

inline std::size_t BinaryBuffer::operator()() const noexcept { return _i; }

inline std::string &BinaryBuffer::reset() noexcept {
  _buffer.clear();
  _i = 0;
  return _buffer;
}

inline void BinaryBuffer::clear() noexcept { std::ignore = reset(); }

inline const std::string &BinaryBuffer::get() const noexcept { return _buffer; }

}  // namespace hictk::hic::internal
