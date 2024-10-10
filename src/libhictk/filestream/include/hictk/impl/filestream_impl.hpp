// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef _WIN32
// clang-format off
#include <windows.h>
#include <errhandlingapi.h>
// clang-format on
#endif

#include <cassert>
#include <cerrno>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iosfwd>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::filestream {

template <typename Mutex>
inline FileStream<Mutex>::FileStream(std::string path, std::shared_ptr<Mutex> mtx,
                                     std::ios::openmode mode)
    : _path(std::move(path)),
      _mtx(std::move(mtx)),
      _ifs(open_file_read(_path, _ifs_flags | std::ios::ate)),
      _ofs(static_cast<bool>(mode & std::ios::out) ? open_file_write(_path, _ofs_flags)
                                                   : std::ofstream{}),
      _file_size(unsafe_tellg()) {
  unsafe_seekg(0);
}

template <typename Mutex>
inline FileStream<Mutex> FileStream<Mutex>::create(std::string path, std::shared_ptr<Mutex> mtx) {
  if (std::filesystem::exists(path)) {
    throw std::runtime_error("file \"" + path + "\" already exists");
  }

  FileStream fs{};
  fs._path = std::move(path);
  fs._mtx = std::move(mtx);
  fs._ofs = open_file_write(fs._path, _ofs_flags | std::ios::trunc);
  fs._ifs = open_file_read(fs._path, _ifs_flags);

  return fs;
}

// We need to explicitly define the move constructor to make things compile on older compilers
// (where std::mutex is not moveable)
template <typename Mutex>
inline FileStream<Mutex>::FileStream(FileStream &&other) noexcept
    : _path(std::move(other._path)),
      _mtx(std::move(other._mtx)),
      _ifs(std::move(other._ifs)),
      _ofs(std::move(other._ofs)),
      _file_size(other._file_size) {}

// We need to explicitly define the move assignment operator to make things compile on older
// compilers (where std::mutex is not moveable)
template <typename Mutex>
inline FileStream<Mutex> &FileStream<Mutex>::operator=(FileStream &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  _path = std::move(other._path);
  _mtx = std::move(other._mtx);
  _ifs = std::move(other._ifs);
  _ofs = std::move(other._ofs);
  _file_size = other._file_size;

  return *this;
}

template <typename Mutex>
inline const std::string &FileStream<Mutex>::path() const noexcept {
  return _path;
}

template <typename Mutex>
inline void FileStream<Mutex>::close() {
  {
    [[maybe_unused]] const auto lck = lock();
    if (_ifs.is_open()) {
      _ifs.close();
    }
    if (_ofs.is_open()) {
      _ofs.close();
    }
  }
  _mtx.reset();
}

template <typename Mutex>
inline void FileStream<Mutex>::seekg(std::streampos position) {
  seekg(position, std::ios::beg);
}

template <typename Mutex>
inline void FileStream<Mutex>::unsafe_seekg(std::streampos position) {
  unsafe_seekg(position, std::ios::beg);
}

template <typename Mutex>
inline void FileStream<Mutex>::seekp(std::streampos position) {
  seekp(position, std::ios::beg);
}

template <typename Mutex>
inline void FileStream<Mutex>::unsafe_seekp(std::streampos position) {
  unsafe_seekp(position, std::ios::beg);
}

template <typename Mutex>
inline void FileStream<Mutex>::seekg(std::streamoff offset, std::ios::seekdir way) {
  [[maybe_unused]] const auto lck = lock();
  unsafe_seekg(offset, way);
}

template <typename Mutex>
inline void FileStream<Mutex>::unsafe_seekg(std::streamoff offset, std::ios::seekdir way) {
  const auto new_pos = new_posg_checked(offset, way);
  _ifs.seekg(new_pos, std::ios::beg);
}

template <typename Mutex>
inline void FileStream<Mutex>::seekp(std::streamoff offset, std::ios::seekdir way) {
  [[maybe_unused]] const auto lck = lock();
  unsafe_seekp(offset, way);
}

template <typename Mutex>
inline void FileStream<Mutex>::unsafe_seekp(std::streamoff offset, std::ios::seekdir way) {
  const auto pos = new_posp(offset, way);
  _ofs.seekp(pos, std::ios::beg);
  _file_size = std::max(static_cast<std::streamsize>(pos), _file_size);
}

template <typename Mutex>
inline std::streampos FileStream<Mutex>::tellg() const {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_tellg();
}

template <typename Mutex>
inline std::streampos FileStream<Mutex>::unsafe_tellg() const {
  const auto offset = _ifs.tellg();
  if (offset >= 0) {
    return offset;
  }
  throw std::runtime_error("FileStream::tellg() called on a bad file handle: " +
                           get_underlying_os_error());
}

template <typename Mutex>
inline std::streampos FileStream<Mutex>::tellp() const {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_tellp();
}

template <typename Mutex>
inline std::streampos FileStream<Mutex>::unsafe_tellp() const {
  const auto offset = _ofs.tellp();
  if (offset >= 0) {
    return offset;
  }
  throw std::runtime_error("FileStream::tellp() called on a bad file handle: " +
                           get_underlying_os_error());
}

template <typename Mutex>
inline std::streamsize FileStream<Mutex>::size() const {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_size();
}

template <typename Mutex>
inline std::streamsize FileStream<Mutex>::unsafe_size() const {
  if constexpr (ndebug_not_defined()) {
    if (path().empty()) {
      return _file_size;
    }
    _ofs.flush();
    const auto offset = unsafe_tellg();
    _ifs.seekg(0, std::ios::end);
    const auto effective_size = unsafe_tellg();
    _ifs.seekg(offset, std::ios::beg);
    if (effective_size != _file_size) {
      throw std::runtime_error("FileStream is corrupted: expected size " +
                               std::to_string(_file_size) + ", found " +
                               std::to_string(effective_size));
    }
  }
  return _file_size;
}

template <typename Mutex>
inline bool FileStream<Mutex>::eof() const {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_eof();
}

template <typename Mutex>
inline bool FileStream<Mutex>::unsafe_eof() const {
  return _ifs.eof();
}

template <typename Mutex>
inline void FileStream<Mutex>::flush() {
  [[maybe_unused]] const auto lck = lock();
  unsafe_flush();
}

template <typename Mutex>
inline void FileStream<Mutex>::unsafe_flush() {
  _ofs.flush();
}

template <typename Mutex>
inline std::unique_lock<Mutex> FileStream<Mutex>::lock() const {
  if (!_mtx) {
    return {};
  }
  return std::unique_lock{*_mtx};
}

template <typename Mutex>
inline bool FileStream<Mutex>::is_locked() const noexcept {
  if (!_mtx) {
    return false;
  }
  // clang8 complains if this expression is written in a single line
  const std::unique_lock lck{*_mtx, std::try_to_lock};
  return !lck.owns_lock();
}

template <typename Mutex>
inline std::string FileStream<Mutex>::read(std::size_t count) {
  std::string buffer(count, '\0');
  read(buffer, count);
  return buffer;
}

template <typename Mutex>
inline void FileStream<Mutex>::read(char *buffer, std::size_t count) {
  [[maybe_unused]] const auto lck = lock();
  unsafe_read(buffer, count);
}

template <typename Mutex>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::seek_and_read(
    std::streamoff offset, char *buffer, std::size_t count, std::ios::seekdir way) {
  auto lck = lock();
  const auto offset1 = unsafe_tellg();
  unsafe_seekg(offset, way);
  unsafe_read(buffer, count);
  // NOLINTNEXTLINE(*-bounds-pointer-arithmetic)
  return std::make_pair(offset1, offset + static_cast<std::streamoff>(count));
}

template <typename Mutex>
inline void FileStream<Mutex>::unsafe_read(char *buffer, std::size_t count) {
  const auto offset1 = unsafe_tellg();
  _ifs.read(buffer, static_cast<std::streamsize>(count));
  const auto bytes_read = unsafe_tellg() - offset1;
  validate_read<char>(bytes_read, count, "FileStream::unsafe_read(char *)");
}

template <typename Mutex>
inline void FileStream<Mutex>::read(std::string &buffer, std::size_t count) {
  try {
    buffer.resize(count);
    if (count == 0) {
      return;
    }
    read(&buffer.front(), count);
  } catch (...) {
    buffer.clear();
    throw;
  }
}

template <typename Mutex>
inline void FileStream<Mutex>::read_append(std::string &buffer, std::size_t count) {
  if (count == 0) {
    return;
  }

  const auto buff_size = buffer.size();
  try {
    buffer.resize(buffer.size() + count);

    [[maybe_unused]] const auto lck = lock();
    const auto offset1 = unsafe_tellg();
    // NOLINTNEXTLINE(*-bounds-pointer-arithmetic)
    _ifs.read(&(*buffer.begin()) + buff_size, static_cast<std::streamsize>(count));
    const auto bytes_read = static_cast<std::size_t>(unsafe_tellg() - offset1);
    validate_read<char>(bytes_read, count, "FileStream::read_append(std::string &)");
  } catch (...) {
    buffer.resize(buff_size);
    throw;
  }
}

template <typename Mutex>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::seek_and_read(
    std::streamoff offset, std::string &buffer, std::size_t count, std::ios::seekdir way) {
  try {
    buffer.resize(count);
    return seek_and_read(offset, buffer.data(), buffer.size(), way);
  } catch (...) {
    buffer.clear();
    throw;
  }
}

template <typename Mutex>
inline bool FileStream<Mutex>::getline(std::string &buffer, char delim) {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_getline(buffer, delim);
}

template <typename Mutex>
inline std::string FileStream<Mutex>::getline(char delim) {
  std::string buffer{};
  getline(buffer, delim);
  return buffer;
}

template <typename Mutex>
inline std::tuple<bool, std::streampos, std::streampos> FileStream<Mutex>::seek_and_getline(
    std::streamoff offset, std::string &buffer, std::ios::seekdir way, char delim) {
  [[maybe_unused]] const auto lck = lock();
  const auto offset1 = unsafe_tellg();
  unsafe_seekg(offset, way);
  const auto status = unsafe_getline(buffer, delim);
  return std::make_tuple(status, offset1, unsafe_tellg());
}

template <typename Mutex>
inline bool FileStream<Mutex>::unsafe_getline(std::string &buffer, char delim) {
  buffer.clear();
  if (unsafe_eof()) {
    _ifs.setstate(std::ios::badbit);
  }
  try {
    return !!std::getline(_ifs, buffer, delim);
  } catch (...) {
    if (unsafe_eof() && !_ifs.bad()) {
      return !!_ifs;
    }
    throw;
  }
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline T FileStream<Mutex>::read() {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_read<T>();
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void FileStream<Mutex>::read(T &buffer) {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_read(buffer);
}

template <typename Mutex>
template <typename Tin, typename Tout, typename std::enable_if_t<std::is_integral_v<Tin>> *>
inline Tout FileStream<Mutex>::read_as_signed() {
  return conditional_static_cast<Tout>(read<Tin>());
}

template <typename Mutex>
template <typename Tin, typename Tout, typename std::enable_if_t<std::is_integral_v<Tin>> *>
inline Tout FileStream<Mutex>::read_as_unsigned() {
  return conditional_static_cast<Tout>(read<Tin>());
}

template <typename Mutex>
template <typename Tin, typename std::enable_if_t<std::is_arithmetic_v<Tin>> *>
inline double FileStream<Mutex>::read_as_double() {
  return conditional_static_cast<double>(read<Tin>());
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void FileStream<Mutex>::unsafe_read(T &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  unsafe_read(reinterpret_cast<char *>(&buffer), sizeof(T));
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline T FileStream<Mutex>::unsafe_read() {
  T buffer{};
  unsafe_read(buffer);
  return buffer;
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void FileStream<Mutex>::read(std::vector<T> &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  read(reinterpret_cast<char *>(&(*buffer.begin())), buffer.size() * sizeof(T));
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline std::vector<T> FileStream<Mutex>::read_vector(std::size_t size) {
  std::vector<T> buffer(size);
  read(buffer);
  return buffer;
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::seek_and_read(
    std::streamoff offset, std::vector<T> &buffer, std::ios::seekdir way) {
  [[maybe_unused]] const auto lck = lock();
  const auto offset1 = unsafe_tellg();
  unsafe_seekg(offset, way);
  unsafe_read(buffer);
  const auto num_bytes = buffer.size() * sizeof(T);
  return std::make_pair(offset1, offset + static_cast<std::streamoff>(num_bytes));
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void FileStream<Mutex>::unsafe_read(std::vector<T> &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  unsafe_read(reinterpret_cast<char *>(&(*buffer.begin())), buffer.size() * sizeof(T));
}

template <typename Mutex>
inline void FileStream<Mutex>::write(const char *buffer, std::size_t count) {
  [[maybe_unused]] const auto lck = lock();
  unsafe_write(buffer, count);
}

template <typename Mutex>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::seek_and_write(
    std::streamoff offset, const char *buffer, std::size_t count, std::ios::seekdir way) {
  [[maybe_unused]] const auto lck = lock();
  const auto offset1 = unsafe_tellp();
  unsafe_seekp(offset, way);
  unsafe_write(buffer, count);
  return std::make_pair(offset1, offset + static_cast<std::streamoff>(count));
}

template <typename Mutex>
inline void FileStream<Mutex>::unsafe_write(const char *buffer, std::size_t count) {
  _ofs.write(buffer, static_cast<std::streamsize>(count));
  _file_size = std::max(static_cast<std::streamsize>(unsafe_tellp()), _file_size);
}

template <typename Mutex>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::append(const char *buffer,
                                                                           std::size_t count) {
  const auto lck = lock();
  unsafe_seekp(0, std::ios::end);
  const auto offset = unsafe_tellp();
  unsafe_write(buffer, count);

  return std::make_pair(offset, offset + static_cast<std::streamoff>(count));
}

template <typename Mutex>
inline void FileStream<Mutex>::write(std::string_view buffer) {
  return write(buffer.data(), buffer.size());
}

template <typename Mutex>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::seek_and_write(
    std::streamoff offset, std::string_view buffer, std::ios::seekdir way) {
  return seek_and_write(offset, buffer.data(), buffer.size(), way);
}

template <typename Mutex>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::append(
    std::string_view buffer) {
  return append(buffer.data(), buffer.size());
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void FileStream<Mutex>::write(T buffer) {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_write(buffer);
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void FileStream<Mutex>::unsafe_write(T buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return unsafe_write(reinterpret_cast<const char *>(&buffer), sizeof(T));
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
std::pair<std::streampos, std::streampos> FileStream<Mutex>::append(T buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return append(reinterpret_cast<const char *>(&buffer), sizeof(T));
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::seek_and_write(
    std::streamoff offset, T buffer, std::ios::seekdir way) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return seek_and_write(offset, reinterpret_cast<const char *>(&buffer), sizeof(T), way);
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline void FileStream<Mutex>::write(const std::vector<T> &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return write(reinterpret_cast<const char *>(buffer.data()), buffer.size() * sizeof(T));
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::seek_and_write(
    std::streamoff offset, const std::vector<T> &buffer, std::ios::seekdir way) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return seek_and_write(offset, reinterpret_cast<const char *>(buffer.data()),
                        buffer.size() * sizeof(T), way);
}

template <typename Mutex>
template <typename T, typename std::enable_if_t<std::is_arithmetic_v<T>> *>
inline std::pair<std::streampos, std::streampos> FileStream<Mutex>::append(
    const std::vector<T> &buffer) {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  return append(reinterpret_cast<const char *>(buffer.data()), buffer.size() * sizeof(T));
}

template <typename Mutex>
inline void FileStream<Mutex>::resize(std::streamsize new_size) {
  [[maybe_unused]] const auto lck = lock();
  if (!_ofs.is_open()) {
    throw std::runtime_error("FileStream::resize() was called on a file opened in read-only mode");
  }

  if (new_size != _file_size) {
    std::filesystem::resize_file(_path, conditional_static_cast<std::uintmax_t>(new_size));
    const auto current_ifs_offset = unsafe_tellg();
    const auto current_ofs_offset = unsafe_tellp();
    _ifs = open_file_read(_path, _ifs_flags);
    _ofs = open_file_write(_path, _ofs_flags);
    _file_size = new_size;

    unsafe_seekg(std::min(conditional_static_cast<std::streampos>(_file_size),
                          conditional_static_cast<std::streampos>(current_ifs_offset)));
    unsafe_seekp(std::min(conditional_static_cast<std::streampos>(_file_size),
                          conditional_static_cast<std::streampos>(current_ofs_offset)));
  }
}

template <typename Mutex>
inline std::streampos FileStream<Mutex>::new_posg(std::streamoff offset, std::ios::seekdir way) {
  switch (way) {
    case std::ios::beg:
      return offset;
    case std::ios::cur:
      return unsafe_tellg() + offset;
    case std::ios::end:
      return static_cast<std::streampos>(_file_size) - offset;
    default:
      HICTK_UNREACHABLE_CODE;
  }
}

template <typename Mutex>
inline std::streampos FileStream<Mutex>::new_posg_checked(std::streamoff offset,
                                                          std::ios::seekdir way) {
  const auto new_pos = new_posg(offset, way);
  if (new_pos < 0 || new_pos >= static_cast<std::streampos>(_file_size + 1)) {
    auto msg1 = "new_posg_checked returned invalid offset=" + std::to_string(new_pos) +
                "; offset not between 0 and " + std::to_string(_file_size + 1);
    if (const auto msg2 = get_underlying_os_error(); msg1 != "Success") {
      msg1 += msg2;
    }

    throw std::runtime_error(msg1);
  }
  return new_pos;
}

template <typename Mutex>
inline std::streampos FileStream<Mutex>::new_posp(std::streamoff offset, std::ios::seekdir way) {
  switch (way) {
    case std::ios::beg:
      return offset;
    case std::ios::cur:
      return unsafe_tellp() + offset;
    case std::ios::end:
      return static_cast<std::streampos>(_file_size) - offset;
    default:
      HICTK_UNREACHABLE_CODE;
  }
}
template <typename Mutex>
inline void FileStream<Mutex>::update_file_size() {
  [[maybe_unused]] const auto lck = lock();
  unsafe_update_file_size();
}

template <typename Mutex>
inline void FileStream<Mutex>::unsafe_update_file_size() {
  unsafe_flush();
  const auto offset = unsafe_tellg();
  unsafe_seekg(0, std::ios::end);
  _file_size = std::max(static_cast<std::streamsize>(unsafe_tellg()), _file_size);
  unsafe_seekg(offset);
}

template <typename Mutex>
inline std::ifstream FileStream<Mutex>::open_file_read(const std::string &path,
                                                       std::ifstream::openmode mode) {
  std::ifstream fs;
  fs.exceptions(fs.exceptions() | std::ios::failbit | std::ios::badbit);
  fs.open(path, mode);
  return fs;
}

template <typename Mutex>
inline std::ofstream FileStream<Mutex>::open_file_write(const std::string &path,
                                                        std::ofstream::openmode mode) {
  std::ofstream fs;
  fs.exceptions(fs.exceptions() | std::ios::failbit | std::ios::badbit);
  fs.open(path, mode);
  return fs;
}

template <typename Mutex>
inline std::string FileStream<Mutex>::get_underlying_os_error() {
  return get_underlying_os_error(errno);
}

template <typename Mutex>
inline std::string FileStream<Mutex>::get_underlying_os_error(int errno_) {
  std::string buffer(256, '\0');  // NOLINT(*-avoid-magic-numbers)
  get_underlying_os_error(errno_, buffer);
  return buffer;
}

template <typename Mutex>
inline void FileStream<Mutex>::get_underlying_os_error([[maybe_unused]] int errno_,
                                                       std::string &buffer) {
#ifdef _WIN32
  if (const auto ec = GetLastError(); ec != 0) {
    LPSTR msg_buffer = nullptr;

    const auto msg_size = FormatMessageA(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        nullptr, ec, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&msg_buffer, 0, nullptr);

    buffer.assign(msg_buffer, msg_size);
    LocalFree(msg_buffer);
    return;
  }
#endif

  if (errno_ == 0) {
    buffer = "Success";
    return;
  }
#if defined(_GNU_SOURCE)
  buffer = strerror_r(errno_, buffer.data(), buffer.size());
#elif defined(_WIN32)
  buffer.resize(std::max(buffer.capacity(), std::size_t{1024}), '\0');
  const int status = strerror_s(buffer.data(), buffer.size(), errno_);
#else
  buffer.resize(std::max(buffer.capacity(), std::size_t{1024}), '\0');
  const int status = strerror_r(errno_, buffer.data(), buffer.size());
#endif

#ifndef _GNU_SOURCE
  switch (status) {
    case 0: {
      // strerror_r/s call was successful
      if (const auto pos = buffer.find('\0'); pos != std::string::npos) {
        if (pos == 0) {
          // this should never happen
          buffer.clear();
        } else {
          buffer.resize(pos);
        }
      }
      return;
    }
    case EINVAL: {
#ifdef _WIN32
      buffer = "strerror_s: unknown errno: " + std::to_string(errno_);
#else
      buffer = "strerror_r: unknown errno: " + std::to_string(errno_);
#endif
      return;
    }
    case ERANGE: {
      buffer.resize(buffer.size() * 2, '\0');
      return get_underlying_os_error(errno_, buffer);
    }
    default: {
      assert(status < 0);
      // on old versions of glibc error is communicated through a new errno value
      get_underlying_os_error(errno, buffer);
      buffer = "strerror_r: " + buffer;
    }
  }
#endif
}

template <typename Mutex>
template <typename T, typename I1, typename I2>
inline void FileStream<Mutex>::validate_read(I1 bytes_read, I2 count_expected,
                                             std::string_view prefix) {
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
  throw std::runtime_error(std::string{prefix} + " failed: expected to read " +
                           std::to_string(bytes_expected) + " bytes, but only read " +
                           std::to_string(bytes_read));
}

}  // namespace hictk::filestream
