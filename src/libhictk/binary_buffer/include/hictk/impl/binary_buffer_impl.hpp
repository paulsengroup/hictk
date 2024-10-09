// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <string>
#include <tuple>
#include <type_traits>

namespace hictk {

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline T BinaryBuffer::read() {
  static_assert(sizeof(char) == 1);
  assert(_i < _buffer.size());
  T x{};

  // NOLINTNEXTLINE(*-pointer-arithmetic,*-type-reinterpret-cast)
  std::memcpy(static_cast<void *>(&x), _buffer.data() + _i, sizeof(T));
  _i += sizeof(T);
  return x;
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void BinaryBuffer::read(T &buff) {
  buff = read<T>();
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void BinaryBuffer::read(std::vector<T> &buff) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  read(reinterpret_cast<char *>(buff.data()), sizeof(T) * buff.size());
}

inline void BinaryBuffer::read(std::string &buff, std::size_t n) {
  buff.resize(n);
  read(buff.data(), n);
}

inline void BinaryBuffer::read(char *buff, std::size_t n) {
  static_assert(sizeof(char) == 1);
  const auto size = n * sizeof(char);
  assert(_i + size < _buffer.size());
  // NOLINTNEXTLINE(*-pointer-arithmetic)
  std::memcpy(static_cast<void *>(buff), _buffer.data() + _i, size);
  _i += size;
}

inline std::string BinaryBuffer::getline(char delim) {
  const std::string_view view{_buffer};
  const auto pos = view.substr(_i).find(delim);
  _i = std::min(pos, _buffer.size());
  return std::string{view.substr(0, pos)};
}

inline void BinaryBuffer::write(const char *data, std::size_t count, bool add_nullterm) {
  _buffer.append(data, count);
  if (add_nullterm) {
    _buffer.append("\0", 1);
  }
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void BinaryBuffer::write(T data) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  write(reinterpret_cast<const char *>(&data), sizeof(T), false);
}

inline void BinaryBuffer::write(const std::string &data, bool add_nullterm) {
  // NOLINTNEXTLINE(*-pointer-arithmetic)
  write(data.c_str(), data.size() + static_cast<std::size_t>(add_nullterm), false);
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void BinaryBuffer::write(const std::vector<T> &data) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  write(reinterpret_cast<const char *>(data.data()), data.size() * sizeof(T), false);
}

inline std::size_t BinaryBuffer::operator()() const noexcept { return _i; }

inline std::string &BinaryBuffer::reset() noexcept {
  _buffer.clear();
  _i = 0;
  return _buffer;
}

inline void BinaryBuffer::clear() noexcept { std::ignore = reset(); }

inline const std::string &BinaryBuffer::get() const noexcept { return _buffer; }

}  // namespace hictk
