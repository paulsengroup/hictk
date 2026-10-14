// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <string>
#include <type_traits>
#include <vector>

namespace hictk {

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
T BinaryBuffer::read() {
  static_assert(sizeof(char) == 1);
  assert(_i < _buffer.size());
  T x{};

  // NOLINTNEXTLINE(*-pointer-arithmetic,*-type-reinterpret-cast)
  std::memcpy(static_cast<void *>(&x), _buffer.data() + _i, sizeof(T));
  _i += sizeof(T);
  return x;
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void BinaryBuffer::read(T &buff) {
  buff = read<T>();
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void BinaryBuffer::read(std::vector<T> &buff) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  read(reinterpret_cast<char *>(buff.data()), sizeof(T) * buff.size());
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>  // NOLINT
void BinaryBuffer::write(T data) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  write(reinterpret_cast<const char *>(&data), sizeof(T), false);
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void BinaryBuffer::write(const std::vector<T> &data) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  write(reinterpret_cast<const char *>(data.data()), data.size() * sizeof(T), false);
}

}  // namespace hictk
