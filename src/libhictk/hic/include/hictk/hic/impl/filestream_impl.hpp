// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hictk/common.hpp"

namespace hictk::hic::internal::filestream {

inline FileStream::FileStream(std::string path)
    : path_(std::move(path)),
      handle_(open_file(path_, std::ios::binary | std::ios::ate)),
      file_size_(static_cast<std::size_t>(handle_.tellg())) {
  handle_.seekg(std::ios::beg);
}

inline const std::string &FileStream::path() const noexcept { return path_; }
inline const std::string &FileStream::url() const noexcept { return path(); }

inline std::size_t FileStream::size() const { return file_size_; }

inline void FileStream::seekg(std::streamoff offset, std::ios::seekdir way) {
  const auto new_pos = this->new_pos(offset, way);
  if (new_pos < 0 || new_pos >= std::int64_t(size() + 1)) {
    throw std::runtime_error("caught an attempt of out-of-bound read");
  }
  handle_.seekg(new_pos, std::ios::beg);
}

inline std::size_t FileStream::tellg() const noexcept {
  return static_cast<std::size_t>(handle_.tellg());
}

inline bool FileStream::eof() const noexcept { return handle_.eof(); }

inline void FileStream::read(std::string &buffer, std::size_t count) {
  buffer.resize(count);
  if (count > 0) {
    return read(&buffer.front(), count);
  }
}

inline void FileStream::read(char *buffer, std::size_t count) {
  handle_.read(buffer, std::int64_t(count));
}

inline void FileStream::append(std::string &buffer, std::size_t count) {
  if (count == 0) {
    return;
  }
  const auto buff_size = buffer.size();
  buffer.resize(buffer.size() + count);

  handle_.read(&(*buffer.begin()) + buff_size, std::int64_t(count));
}

inline bool FileStream::getline(std::string &buffer, char delim) {
  buffer.clear();
  if (eof()) {
    handle_.setstate(std::ios::badbit);
  }
  try {
    return !!std::getline(handle_, buffer, delim);
  } catch (const std::exception &) {
    if (handle_.eof() && !handle_.bad()) {
      return !!handle_;
    }
    throw;
  }
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void FileStream::read(T &buffer) {
  static_assert(sizeof(char) == 1, "");
  return read(reinterpret_cast<char *>(&buffer), sizeof(T));
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T FileStream::read() {
  T buffer{};
  read(buffer);
  return buffer;
}

template <typename Tin, typename Tout,
          typename std::enable_if<std::is_integral<Tin>::value>::type *>
inline Tout FileStream::read_as_signed() {
  auto tmp = read<Tin>();
  return conditional_static_cast<Tout>(read<Tin>());
}

template <typename Tin, typename Tout,
          typename std::enable_if<std::is_integral<Tin>::value>::type *>
inline Tout FileStream::read_as_unsigned() {
  return conditional_static_cast<Tout>(read<Tin>());
}

template <typename Tin, typename std::enable_if<std::is_arithmetic<Tin>::value>::type *>
inline double FileStream::read_as_double() {
  return conditional_static_cast<double>(read<Tin>());
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void FileStream::read(std::vector<T> &buffer) {
  static_assert(sizeof(char) == 1, "");
  return read(reinterpret_cast<char *>(&(*buffer.begin())), buffer.size() * sizeof(T));
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline std::vector<T> FileStream::read(std::size_t size) {
  assert(size != 0);
  std::vector<T> buffer(size);
  read(buffer);
  return buffer;
}

inline std::string FileStream::getline(char delim) {
  std::string buffer{};
  getline(buffer, delim);
  return buffer;
}

inline std::streampos FileStream::new_pos(std::streamoff offset, std::ios::seekdir way) {
  switch (way) {
    case std::ios::beg:
      return static_cast<std::streampos>(offset);
    case std::ios::cur:
      return std::int64_t(tellg()) + offset;
    case std::ios::end:
      return std::int64_t(file_size_) + offset;
    default:
      HICTK_UNREACHABLE_CODE;
  }
}

inline std::ifstream FileStream::open_file(const std::string &path, std::ifstream::openmode mode) {
  std::ifstream ifs;
  ifs.exceptions(ifs.exceptions() | std::ios::failbit | std::ios::badbit);
  ifs.open(path, mode);
  return ifs;
}

}  // namespace hictk::hic::internal::filestream
