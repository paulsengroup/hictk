// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/filestream.hpp"

#ifdef _WIN32
// clang-format off
#include <windows.h>
#include <errhandlingapi.h>
// clang-format on
#endif

#include <fmt/format.h>

#include <algorithm>
#include <cerrno>
#include <cstddef>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

#include "hictk/common.hpp"

namespace hictk::filestream {

static void shrink_to_fit_error_buffer(std::string &buffer) {
  const auto pos = buffer.find('\0');

  if (pos == std::string::npos) {
    return;
  }
  if (pos > 0) {
    buffer.resize(pos);
    return;
  }

  // this should never happen
  buffer.clear();
}

#ifdef _GNU_SOURCE

// https://stackoverflow.com/questions/41953104/strerror-r-is-incorrectly-declared-on-alpine-linux
[[maybe_unused]] static constexpr char *gnu_strerror_r_wrapper_helper(
    [[maybe_unused]] int result, char *buffer, [[maybe_unused]] int err) noexcept {
  return buffer;
}

[[maybe_unused]] static constexpr char *gnu_strerror_r_wrapper_helper(
    char *result, [[maybe_unused]] char *buffer, [[maybe_unused]] int err) noexcept {
  return result;
}

[[maybe_unused]] static void gnu_strerror_r_wrapper(int errnum, std::string &buff) {
  const char *msg = gnu_strerror_r_wrapper_helper(strerror_r(errnum, buff.data(), buff.size()),
                                                  buff.data(), errnum);
  if (buff.data() == msg) {
    shrink_to_fit_error_buffer(buff);
  } else {
    buff.assign(msg);
  }
}
#endif

FileStream::FileStream(std::string path, std::shared_ptr<Mutex> mtx, std::ios::openmode mode)
    : _path(std::move(path)),
      _mtx(std::move(mtx)),
      _ifs(open_file_read(_path, _ifs_flags | std::ios::ate)),
      _ofs(static_cast<bool>(mode & std::ios::out) ? open_file_write(_path, _ofs_flags)
                                                   : std::ofstream{}),
      _file_size(unsafe_tellg()) {
  unsafe_seekg(0);
}

FileStream FileStream::create(std::string path, std::shared_ptr<Mutex> mtx) {
  if (std::filesystem::exists(path)) {
    throw std::runtime_error(fmt::format(FMT_STRING("file \"{}\" already exists"), path));
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

FileStream::FileStream(FileStream &&other) noexcept
    : _path(std::move(other._path)),
      _mtx(std::move(other._mtx)),
      _ifs(std::move(other._ifs)),
      _ofs(std::move(other._ofs)),
      _file_size(other._file_size) {}

// We need to explicitly define the move assignment operator to make things compile on older
// compilers (where std::mutex is not moveable)

FileStream &FileStream::operator=(FileStream &&other) noexcept {
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

const std::string &FileStream::path() const noexcept { return _path; }

void FileStream::close() {
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

void FileStream::seekg(std::streampos position) { seekg(position, std::ios::beg); }

void FileStream::unsafe_seekg(std::streampos position) { unsafe_seekg(position, std::ios::beg); }

void FileStream::seekp(std::streampos position) { seekp(position, std::ios::beg); }

void FileStream::unsafe_seekp(std::streampos position) { unsafe_seekp(position, std::ios::beg); }

void FileStream::seekg(std::streamoff offset, std::ios::seekdir way) {
  [[maybe_unused]] const auto lck = lock();
  unsafe_seekg(offset, way);
}

void FileStream::unsafe_seekg(std::streamoff offset, std::ios::seekdir way) {
  const auto new_pos = new_posg_checked(offset, way);
  _ifs.seekg(new_pos, std::ios::beg);
}

void FileStream::seekp(std::streamoff offset, std::ios::seekdir way) {
  [[maybe_unused]] const auto lck = lock();
  unsafe_seekp(offset, way);
}

void FileStream::unsafe_seekp(std::streamoff offset, std::ios::seekdir way) {
  const auto pos = new_posp(offset, way);
  _ofs.seekp(pos, std::ios::beg);
  _file_size = std::max(static_cast<std::streamsize>(pos), _file_size);
}

std::streampos FileStream::tellg() const {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_tellg();
}

std::streampos FileStream::unsafe_tellg() const {
  const auto offset = _ifs.tellg();
  if (offset >= 0) {
    return offset;
  }
  throw std::runtime_error(
      fmt::format(FMT_STRING("FileStream::tellg() called on a bad file handle: {}"),
                  get_underlying_os_error()));
}

std::streampos FileStream::tellp() const {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_tellp();
}

std::streampos FileStream::unsafe_tellp() const {
  const auto offset = _ofs.tellp();
  if (offset >= 0) {
    return offset;
  }
  throw std::runtime_error(
      fmt::format(FMT_STRING("FileStream::tellp() called on a bad file handle: {}"),
                  get_underlying_os_error()));
}

std::streamsize FileStream::size() const {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_size();
}

std::streamsize FileStream::unsafe_size() const {
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
      throw std::runtime_error(fmt::format(
          FMT_STRING("FileStream::size(): file \"{}\" is corrupted: expected size {}, found {}"),
          _path, _file_size, static_cast<std::int64_t>(effective_size)));
    }
  }
  return _file_size;
}

bool FileStream::eof() const {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_eof();
}

bool FileStream::unsafe_eof() const { return _ifs.eof(); }

void FileStream::flush() {
  [[maybe_unused]] const auto lck = lock();
  unsafe_flush();
}

void FileStream::unsafe_flush() { _ofs.flush(); }

auto FileStream::lock() const -> std::unique_lock<Mutex> {
  if (!_mtx) {
    return {};
  }
  return std::unique_lock{*_mtx};
}

bool FileStream::is_locked() const noexcept {
  if (!_mtx) {
    return false;
  }
  // clang8 complains if this expression is written in a single line
  const std::unique_lock lck{*_mtx, std::try_to_lock};
  return !lck.owns_lock();
}

std::string FileStream::read(std::size_t count) {
  std::string buffer(count, '\0');
  read(buffer, count);
  return buffer;
}

void FileStream::read(char *buffer, std::size_t count) {
  [[maybe_unused]] const auto lck = lock();
  unsafe_read(buffer, count);
}

std::pair<std::streampos, std::streampos> FileStream::seek_and_read(std::streamoff offset,
                                                                    char *buffer, std::size_t count,
                                                                    std::ios::seekdir way) {
  auto lck = lock();
  const auto offset1 = unsafe_tellg();
  unsafe_seekg(offset, way);
  unsafe_read(buffer, count);
  // NOLINTNEXTLINE(*-bounds-pointer-arithmetic)
  return std::make_pair(offset1, offset + static_cast<std::streamoff>(count));
}

void FileStream::unsafe_read(char *buffer, std::size_t count) {
  const auto offset1 = unsafe_tellg();
  _ifs.read(buffer, static_cast<std::streamsize>(count));
  const auto bytes_read = unsafe_tellg() - offset1;
  validate_read<char>(bytes_read, count, "FileStream::unsafe_read(char *)");
}

void FileStream::read(std::string &buffer, std::size_t count) {
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

void FileStream::read_append(std::string &buffer, std::size_t count) {
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

std::pair<std::streampos, std::streampos> FileStream::seek_and_read(std::streamoff offset,
                                                                    std::string &buffer,
                                                                    std::size_t count,
                                                                    std::ios::seekdir way) {
  try {
    buffer.resize(count);
    return seek_and_read(offset, buffer.data(), buffer.size(), way);
  } catch (...) {
    buffer.clear();
    throw;
  }
}

bool FileStream::getline(std::string &buffer, char delim) {
  [[maybe_unused]] const auto lck = lock();
  return unsafe_getline(buffer, delim);
}

std::string FileStream::getline(char delim) {
  std::string buffer{};
  getline(buffer, delim);
  return buffer;
}

std::tuple<bool, std::streampos, std::streampos> FileStream::seek_and_getline(std::streamoff offset,
                                                                              std::string &buffer,
                                                                              std::ios::seekdir way,
                                                                              char delim) {
  [[maybe_unused]] const auto lck = lock();
  const auto offset1 = unsafe_tellg();
  unsafe_seekg(offset, way);
  const auto status = unsafe_getline(buffer, delim);
  return std::make_tuple(status, offset1, unsafe_tellg());
}

bool FileStream::unsafe_getline(std::string &buffer, char delim) {
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

void FileStream::write(const char *buffer, std::size_t count) {
  [[maybe_unused]] const auto lck = lock();
  unsafe_write(buffer, count);
}

std::pair<std::streampos, std::streampos> FileStream::seek_and_write(std::streamoff offset,
                                                                     const char *buffer,
                                                                     std::size_t count,
                                                                     std::ios::seekdir way) {
  [[maybe_unused]] const auto lck = lock();
  const auto offset1 = unsafe_tellp();
  unsafe_seekp(offset, way);
  unsafe_write(buffer, count);
  return std::make_pair(offset1, offset + static_cast<std::streamoff>(count));
}

void FileStream::unsafe_write(const char *buffer, std::size_t count) {
  _ofs.write(buffer, static_cast<std::streamsize>(count));
  _file_size = std::max(static_cast<std::streamsize>(unsafe_tellp()), _file_size);
}

std::pair<std::streampos, std::streampos> FileStream::append(const char *buffer,
                                                             std::size_t count) {
  const auto lck = lock();
  unsafe_seekp(0, std::ios::end);
  const auto offset = unsafe_tellp();
  unsafe_write(buffer, count);

  return std::make_pair(offset, offset + static_cast<std::streamoff>(count));
}

void FileStream::write(std::string_view buffer) { write(buffer.data(), buffer.size()); }

std::pair<std::streampos, std::streampos> FileStream::seek_and_write(std::streamoff offset,
                                                                     std::string_view buffer,
                                                                     std::ios::seekdir way) {
  return seek_and_write(offset, buffer.data(), buffer.size(), way);
}

std::pair<std::streampos, std::streampos> FileStream::append(std::string_view buffer) {
  return append(buffer.data(), buffer.size());
}

void FileStream::resize(std::streamsize new_size) {
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

std::streampos FileStream::new_posg(std::streamoff offset, std::ios::seekdir way) const {
  switch (way) {
    case std::ios::beg:
      return offset;
    case std::ios::cur:
      return unsafe_tellg() + offset;
    case std::ios::end:
      return static_cast<std::streampos>(_file_size) - offset;
    default:
      unreachable_code();
  }
}

std::streampos FileStream::new_posg_checked(std::streamoff offset, std::ios::seekdir way) const {
  const auto new_pos = new_posg(offset, way);
  if (new_pos < 0 || new_pos >= static_cast<std::streampos>(_file_size + 1)) {
    const auto error = get_underlying_os_error();
    throw std::runtime_error(fmt::format(
        FMT_STRING("new_posg_checked returned invalid offset={}; offset not between 0 and {}{}{}"),
        static_cast<std::int64_t>(new_pos), conditional_static_cast<std::int64_t>(_file_size + 1),
        error.empty() ? "" : ": ", error));
  }
  return new_pos;
}

std::streampos FileStream::new_posp(std::streamoff offset, std::ios::seekdir way) const {
  switch (way) {
    case std::ios::beg:
      return offset;
    case std::ios::cur:
      return unsafe_tellp() + offset;
    case std::ios::end:
      return static_cast<std::streampos>(_file_size) - offset;
    default:
      unreachable_code();
  }
}

std::ifstream FileStream::open_file_read(const std::string &path, std::ifstream::openmode mode) {
  std::ifstream fs;
  fs.exceptions(fs.exceptions() | std::ios::failbit | std::ios::badbit);
  fs.open(path, mode);
  return fs;
}

std::ofstream FileStream::open_file_write(const std::string &path, std::ofstream::openmode mode) {
  std::ofstream fs;
  fs.exceptions(fs.exceptions() | std::ios::failbit | std::ios::badbit);
  fs.open(path, mode);
  return fs;
}

std::string FileStream::get_underlying_os_error() { return get_underlying_os_error(errno); }

std::string FileStream::get_underlying_os_error(int errno_) {
  std::string buffer;
  get_underlying_os_error(errno_, buffer);
  return buffer;
}

static void initialize_error_buffer(std::string &s, std::size_t min_size = 1024) {
#ifdef _WIN32
  // This buffer cannot be larger than 64K bytes
  // https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-formatmessagea
  min_size = std::min(min_size, 64 * 1024ULL);
#endif
  s.clear();
  s.resize(std::max(min_size, s.capacity()), '\0');
}

// NOLINTBEGIN(*-use-concise-preprocessor-directives)
#if defined(_WIN32)
[[nodiscard]] static bool get_last_error_helper(std::string &buffer) noexcept {
  if (const auto ec = GetLastError(); ec != 0) {
    auto msg_size = FormatMessageA(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
                                   nullptr, ec, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
                                   conditional_static_cast<LPSTR>(buffer.data()),
                                   conditional_static_cast<DWORD>(buffer.size()), nullptr);

    if (msg_size == 0) {
      return false;
    }

    // for some reason, certain messages end with a newline...
    if (msg_size > 2 && buffer[msg_size - 2] == '\r' && buffer[msg_size - 1] == '\n') {
      msg_size -= 2;
    }

    assert(msg_size <= buffer.capacity());
    buffer.resize(msg_size);
    return true;
  }
  return false;
}

[[nodiscard]] static int get_underlying_os_error_helper(int errno_, std::string &buffer) noexcept {
  if (get_last_error_helper(buffer)) {
    return 0;
  }
  return strerror_s(buffer.data(), buffer.size(), errno_);
}
#elif defined(_GNU_SOURCE)
[[nodiscard]] static int get_underlying_os_error_helper(int errno_, std::string &buffer) noexcept {
  gnu_strerror_r_wrapper(errno_, buffer);
  return 0;
}
#else
[[nodiscard]] static int get_underlying_os_error_helper(int errno_, std::string &buffer) noexcept {
  return strerror_r(errno_, buffer.data(), buffer.size());
}
#endif
// NOLINTEND(*-use-concise-preprocessor-directives)

void FileStream::get_underlying_os_error(int errno_, std::string &buffer) {
  initialize_error_buffer(buffer);
  const auto status = get_underlying_os_error_helper(errno_, buffer);

  switch (status) {
    case 0: {
      // strerror_r/s call was successful
      shrink_to_fit_error_buffer(buffer);
      return;
    }
    case EINVAL: {
      buffer = fmt::format(FMT_STRING("{}: unknown errno: {}"),
#ifdef _WIN32
                           "strerror_s",
#else
                           "strerror_r",
#endif
                           errno_);
      return;
    }
    case ERANGE: {
      initialize_error_buffer(buffer, buffer.capacity() * 2);
      get_underlying_os_error(errno_, buffer);
      return;
    }
    default: {
      assert(status < 0);
      // on old versions of glibc error is communicated through a new errno value
      get_underlying_os_error(errno, buffer);
      buffer = fmt::format(FMT_STRING("strerror_r failed: {}"), buffer);
    }
  }
}

void FileStream::raise_read_validation_error(std::string_view prefix, std::int64_t bytes_expected,
                                             std::int64_t bytes_read) {
  throw std::runtime_error(
      fmt::format(FMT_STRING("{} failed: expected to read {} bytes, but only read {}"), prefix,
                  bytes_expected, bytes_read));
}

}  // namespace hictk::filestream
