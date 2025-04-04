// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <archive.h>
#include <fmt/format.h>
#include <fmt/std.h>

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>

#include "hictk/common.hpp"
#include "hictk/tools/compressed_io.hpp"

namespace hictk::tools::io {

CompressedReader::CompressedReader(const std::filesystem::path& path, std::size_t buff_capacity) {
  _buff.reserve(buff_capacity);
  open(path);
}

void CompressedReader::open(const std::filesystem::path& path) {
  auto handle_open_errors = [&](la_ssize_t status) {
    if (status == ARCHIVE_EOF) {
      _eof = true;
    }
    if (status < ARCHIVE_OK) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to open file {} for reading (error code {}): {}"), _path,
                      archive_errno(_arc.get()), archive_error_string(_arc.get())));
    }
  };

  if (is_open()) {
    close();
  }

  _path = path;
  if (_path.empty()) {
    return;
  }

  _arc.reset(archive_read_new());
  if (!_arc) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to allocate a buffer of to read file {}"), _path));
  }

  handle_open_errors(archive_read_support_filter_all(_arc.get()));
  handle_open_errors(archive_read_support_format_empty(_arc.get()));
  handle_open_errors(archive_read_support_format_raw(_arc.get()));
#ifdef _MSC_VER
  handle_open_errors(archive_read_open_filename_w(_arc.get(), _path.c_str(), _buff.capacity()));
#else
  handle_open_errors(archive_read_open_filename(_arc.get(), _path.c_str(), _buff.capacity()));
#endif
  handle_open_errors(archive_read_next_header(_arc.get(), _arc_entry.get()));
  _idx = 0;
}

CompressedReader::operator bool() const { return is_open() && !eof(); }

bool CompressedReader::operator!() const { return !operator bool(); }

bool CompressedReader::eof() const noexcept {
  assert(is_open());
  return _eof;
}

bool CompressedReader::is_open() const noexcept { return !!_arc; }

void CompressedReader::close() {
  if (is_open()) {
    _arc = nullptr;
    _buff.clear();
    _eof = false;
  }
}

void CompressedReader::reset() {
  close();
  open(_path);
  _idx = 0;
  _buff.clear();
  _tok_tmp_buff.clear();
}

const std::filesystem::path& CompressedReader::path() const noexcept { return _path; }
std::string CompressedReader::path_string() const noexcept { return _path.string(); }

void CompressedReader::handle_libarchive_errors(la_ssize_t errcode) const {
  if (errcode < ARCHIVE_OK) {
    handle_libarchive_errors();
  }
}

void CompressedReader::handle_libarchive_errors() const {
  if (const auto status = archive_errno(_arc.get()); status < ARCHIVE_OK) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("the following error occurred while reading file (error code {}): {}"), _path,
        archive_errno(_arc.get()), archive_error_string(_arc.get())));
  }
}

bool CompressedReader::getline(std::string& buff, char sep) {
  assert(is_open());
  buff.clear();
  if (eof()) {
    return false;
  }

  while (!read_next_token(buff, sep)) {
    if (!read_next_chunk()) {
      assert(eof());
      return !buff.empty();
    }
  }
  return true;
}

std::string_view CompressedReader::getline(char sep) {
  assert(is_open());
  if (eof()) {
    return std::string_view{};
  }

  _tok_tmp_buff.clear();
  while (true) {
    if (const auto tok = read_next_token(sep); !tok.empty()) {
      return tok;  // NOLINT
    }
    if (!read_next_chunk()) {
      assert(eof());
      return std::string_view{};
    }
  }
}

bool CompressedReader::readall(std::string& buff, char sep) {
  assert(is_open());
  buff.clear();
  if (eof()) {
    return false;
  }

  while (!eof()) {  // NOLINT
    std::ignore = getline(sep);
    buff.append(_buff.begin(), _buff.end());
    _idx = 0;
    _buff.clear();
  }
  return true;
}

std::string CompressedReader::readall(char sep) {
  std::string buff;
  readall(buff, sep);
  return buff;
}

bool CompressedReader::read_next_chunk() {
  assert(!eof());
  assert(is_open());
  _buff.resize(_buff.capacity());
  const auto bytes_read = archive_read_data(_arc.get(), _buff.data(), _buff.capacity());
  if (bytes_read < 0) {
    handle_libarchive_errors();
  } else if (bytes_read == 0) {
    _eof = true;
    _buff.clear();
    _tok_tmp_buff.clear();
    return false;
  }
  _buff.resize(static_cast<std::size_t>(bytes_read));
  _idx = 0;
  return true;
}

bool CompressedReader::read_next_token(std::string& buff, char sep) {
  assert(!eof());
  assert(is_open());
  assert(_idx <= _buff.size());
  if (_idx == _buff.size()) {
    return false;
  }

  const auto pos = _buff.find(sep, _idx);
  const auto i = static_cast<std::int64_t>(_idx);
  if (pos == std::string::npos) {
    buff.append(_buff.begin() + i, _buff.end());
    return false;
  }

  assert(pos >= _idx);
  buff.append(_buff.begin() + i, _buff.begin() + static_cast<std::int64_t>(pos));
  _idx = pos + 1;
  return true;
}

std::string_view CompressedReader::read_next_token(char sep) {
  assert(!eof());
  assert(is_open());
  assert(_idx <= _buff.size());
  if (_idx == _buff.size()) {
    return std::string_view{};
  }

  const auto pos = _buff.find(sep, _idx);
  const auto i = static_cast<std::int64_t>(_idx);
  if (pos == std::string::npos) {
    _tok_tmp_buff.append(_buff.begin() + i, _buff.end());
    return std::string_view{};
  }

  assert(pos >= _idx);
  _idx = pos + 1;
  if (_tok_tmp_buff.empty()) {  // NOLINTNEXTLINE
    return std::string_view{_buff.data() + static_cast<std::size_t>(i),
                            pos - static_cast<std::size_t>(i)};
  }

  _tok_tmp_buff.append(_buff.begin() + i, _buff.begin() + static_cast<std::int64_t>(pos));
  return std::string_view{_tok_tmp_buff};
}
}  // namespace hictk::tools::io
