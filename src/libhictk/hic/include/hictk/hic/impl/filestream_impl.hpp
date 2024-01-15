// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iosfwd>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/common.hpp"

namespace hictk::hic::internal::filestream {

inline FileStream::FileStream(std::string path)
    : _path(std::move(path)),
      _ifs(open_file_read(_path, std::ios::binary | std::ios::ate)),
      _file_size(static_cast<std::size_t>(_ifs.tellg())) {
  _ifs.seekg(0, std::ios::beg);
}

inline FileStream FileStream::create(std::string path) {
  if (std::filesystem::exists(path)) {
    throw std::runtime_error("file\"" + path + "\" already exists");
  }
  FileStream fs{};
  fs._path = std::move(path);
  fs._ofs = open_file_write(fs._path, std::ios::trunc | std::ios::binary);
  fs._ifs = open_file_read(fs._path, std::ios::binary);

  return fs;
}

inline const std::string &FileStream::path() const noexcept { return _path; }
inline const std::string &FileStream::url() const noexcept { return path(); }

inline std::size_t FileStream::size() const { return _file_size; }

inline void FileStream::seekg(std::streamoff offset, std::ios::seekdir way) {
  const auto new_pos = new_posg(offset, way);
  if (new_pos < 0 || new_pos >= std::int64_t(size() + 1)) {
    throw std::runtime_error("caught an attempt of out-of-bound read");
  }
  _ifs.seekg(new_pos, std::ios::beg);
}

inline std::size_t FileStream::tellg() const noexcept {
  return static_cast<std::size_t>(_ifs.tellg());
}

inline void FileStream::seekp(std::streamoff offset, std::ios::seekdir way) {
  _ofs.seekp(new_posp(offset, way), std::ios::beg);
}

inline std::size_t FileStream::tellp() const noexcept {
  return static_cast<std::size_t>(_ofs.tellp());
}

inline bool FileStream::eof() const noexcept { return _ifs.eof(); }

inline void FileStream::flush() { _ofs.flush(); }

inline void FileStream::read(std::string &buffer, std::size_t count) {
  buffer.resize(count);
  if (count > 0) {
    return read(&buffer.front(), count);
  }
}

inline void FileStream::read(char *buffer, std::size_t count) {
  _ifs.read(buffer, std::int64_t(count));
}

inline void FileStream::read_append(std::string &buffer, std::size_t count) {
  if (count == 0) {
    return;
  }
  const auto buff_size = buffer.size();
  buffer.resize(buffer.size() + count);

  _ifs.read(&(*buffer.begin()) + buff_size, std::int64_t(count));
}

inline bool FileStream::getline(std::string &buffer, char delim) {
  buffer.clear();
  if (eof()) {
    _ifs.setstate(std::ios::badbit);
  }
  try {
    return !!std::getline(_ifs, buffer, delim);
  } catch (const std::exception &) {
    if (_ifs.eof() && !_ifs.bad()) {
      return !!_ifs;
    }
    throw;
  }
}

inline void FileStream::write(std::string_view buffer) {
  return write(buffer.data(), buffer.size());
}

inline void FileStream::write(const char *buffer, std::size_t count) {
  _ofs.write(buffer, std::int64_t(count));
  _file_size = std::max(static_cast<std::size_t>(_ofs.tellp()), _file_size);
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void FileStream::read(T &buffer) {
  static_assert(sizeof(char) == 1);
  return read(reinterpret_cast<char *>(&buffer), sizeof(T));
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T FileStream::read() {
  T buffer{};
  read(buffer);
  return buffer;
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void FileStream::write(T buffer) {
  static_assert(sizeof(char) == 1);
  return write(reinterpret_cast<const char *>(&buffer), sizeof(T));
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
  static_assert(sizeof(char) == 1);
  return read(reinterpret_cast<char *>(&(*buffer.begin())), buffer.size() * sizeof(T));
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline void FileStream::write(const std::vector<T> &buffer) {
  static_assert(sizeof(char) == 1);
  return write(reinterpret_cast<const char *>(buffer.data()), buffer.size() * sizeof(T));
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

inline std::streampos FileStream::new_posg(std::streamoff offset, std::ios::seekdir way) {
  switch (way) {
    case std::ios::beg:
      return static_cast<std::streampos>(offset);
    case std::ios::cur:
      return std::int64_t(tellg()) + offset;
    case std::ios::end:
      return std::int64_t(_file_size) + offset;
    default:
      HICTK_UNREACHABLE_CODE;
  }
}

inline std::streampos FileStream::new_posp(std::streamoff offset, std::ios::seekdir way) {
  switch (way) {
    case std::ios::beg:
      return static_cast<std::streampos>(offset);
    case std::ios::cur:
      return std::int64_t(tellp()) + offset;
    case std::ios::end:
      return std::int64_t(_file_size) + offset;
    default:
      HICTK_UNREACHABLE_CODE;
  }
}

inline void FileStream::update_file_size() {
  const auto offset = _ifs.tellg();
  _ifs.seekg(0, std::ios::end);
  _file_size = std::max(static_cast<std::size_t>(_ifs.tellg()), _file_size);
  _ifs.seekg(offset, std::ios::beg);
}

inline std::ifstream FileStream::open_file_read(const std::string &path,
                                                std::ifstream::openmode mode) {
  std::ifstream fs;
  fs.exceptions(fs.exceptions() | std::ios::failbit | std::ios::badbit);
  fs.open(path, mode);
  return fs;
}

inline std::ofstream FileStream::open_file_write(const std::string &path,
                                                 std::ofstream::openmode mode) {
  std::ofstream fs;
  fs.exceptions(fs.exceptions() | std::ios::failbit | std::ios::badbit);
  fs.open(path, mode);
  return fs;
}

}  // namespace hictk::hic::internal::filestream
