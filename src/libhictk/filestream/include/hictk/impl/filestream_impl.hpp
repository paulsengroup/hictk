// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <ios>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::filestream {

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
T FileStream::read() {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_read<T>();
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void FileStream::read(T &buffer) {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_read(buffer);
}

template <typename Tin, typename Tout, typename std::enable_if_t<std::is_integral_v<Tin>> *>
Tout FileStream::read_as_signed() {
  return conditional_static_cast<Tout>(read<Tin>());
}

template <typename Tin, typename Tout, typename std::enable_if_t<std::is_integral_v<Tin>> *>
Tout FileStream::read_as_unsigned() {
  return conditional_static_cast<Tout>(read<Tin>());
}

template <typename Tin, typename std::enable_if_t<std::is_arithmetic_v<Tin>> *>
double FileStream::read_as_double() {
  return conditional_static_cast<double>(read<Tin>());
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void FileStream::unsafe_read(T &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  unsafe_read(reinterpret_cast<char *>(&buffer), sizeof(T));
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
T FileStream::unsafe_read() {
  T buffer{};
  unsafe_read(buffer);
  return buffer;
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void FileStream::read(std::vector<T> &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  read(reinterpret_cast<char *>(&(*buffer.begin())), buffer.size() * sizeof(T));
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
std::vector<T> FileStream::read_vector(std::size_t size) {
  std::vector<T> buffer(size);
  read(buffer);
  return buffer;
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
std::pair<std::streampos, std::streampos> FileStream::seek_and_read(std::streamoff offset,
                                                                    std::vector<T> &buffer,
                                                                    std::ios::seekdir way) {
  [[maybe_unused]] const auto lck = lock();
  const auto offset1 = unsafe_tellg();
  unsafe_seekg(offset, way);
  unsafe_read(buffer);
  const auto num_bytes = buffer.size() * sizeof(T);
  return std::make_pair(offset1, offset + static_cast<std::streamoff>(num_bytes));
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void FileStream::unsafe_read(std::vector<T> &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  unsafe_read(reinterpret_cast<char *>(&(*buffer.begin())), buffer.size() * sizeof(T));
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void FileStream::write(T buffer) {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_write(buffer);
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void FileStream::unsafe_write(T buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return unsafe_write(reinterpret_cast<const char *>(&buffer), sizeof(T));
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
std::pair<std::streampos, std::streampos> FileStream::append(T buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return append(reinterpret_cast<const char *>(&buffer), sizeof(T));
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
std::pair<std::streampos, std::streampos> FileStream::seek_and_write(std::streamoff offset,
                                                                     T buffer,
                                                                     std::ios::seekdir way) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return seek_and_write(offset, reinterpret_cast<const char *>(&buffer), sizeof(T), way);
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
void FileStream::write(const std::vector<T> &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return write(reinterpret_cast<const char *>(buffer.data()), buffer.size() * sizeof(T));
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
std::pair<std::streampos, std::streampos> FileStream::seek_and_write(std::streamoff offset,
                                                                     const std::vector<T> &buffer,
                                                                     std::ios::seekdir way) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return seek_and_write(offset, reinterpret_cast<const char *>(buffer.data()),
                        buffer.size() * sizeof(T), way);
}

template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
std::pair<std::streampos, std::streampos> FileStream::append(const std::vector<T> &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return append(reinterpret_cast<const char *>(buffer.data()), buffer.size() * sizeof(T));
}

template <typename T, typename I1, typename I2>
void FileStream::validate_read(I1 bytes_read, I2 count_expected, std::string_view prefix) {
  static_assert(std::is_arithmetic_v<T>);
  static_assert(std::is_integral_v<I1> || is_specialization_v<I1, std::fpos>);
  static_assert(std::is_integral_v<I2> || is_specialization_v<I1, std::fpos>);

  if constexpr (std::is_signed_v<I1>) {
    assert(bytes_read >= 0);
  }
  if constexpr (std::is_signed_v<I2>) {
    assert(count_expected >= 0);
  }

  const auto bytes_expected = conditional_static_cast<std::size_t>(count_expected) * sizeof(T);
  if (conditional_static_cast<std::size_t>(bytes_read) == bytes_expected) {
    return;
  }

  assert(conditional_static_cast<std::size_t>(bytes_read) < bytes_expected);
  raise_read_validation_error(prefix, bytes_expected, bytes_read);
}

}  // namespace hictk::filestream
