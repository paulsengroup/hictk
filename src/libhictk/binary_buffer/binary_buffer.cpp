// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/binary_buffer.hpp"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>

namespace hictk {

void BinaryBuffer::read(std::string &buff, std::size_t n) {
  buff.resize(n);
  read(buff.data(), n);
}

void BinaryBuffer::read(char *buff, std::size_t n) {
  static_assert(sizeof(char) == 1);
  const auto size = n * sizeof(char);
  assert(_i + size <= _buffer.size());
  // NOLINTNEXTLINE(*-pointer-arithmetic)
  std::memcpy(static_cast<void *>(buff), _buffer.data() + _i, size);
  _i += size;
}

std::string BinaryBuffer::getline(char delim) {
  const std::string_view view{_buffer};
  const auto pos = view.substr(_i).find(delim);
  _i = std::min(pos, _buffer.size());
  return std::string{view.substr(0, pos)};
}

void BinaryBuffer::write(const char *data, std::size_t count, bool add_nullterm) {
  _buffer.append(data, count);
  if (add_nullterm) {
    _buffer.append("\0", 1);
  }
}

void BinaryBuffer::write(const std::string &data, bool add_nullterm) {
  // NOLINTNEXTLINE(*-pointer-arithmetic)
  write(data.c_str(), data.size() + static_cast<std::size_t>(add_nullterm), false);
}

std::size_t BinaryBuffer::operator()() const noexcept { return _i; }

std::string &BinaryBuffer::reset() noexcept {
  _buffer.clear();
  _i = 0;
  return _buffer;
}

void BinaryBuffer::clear() noexcept { std::ignore = reset(); }

const std::string &BinaryBuffer::get() const noexcept { return _buffer; }

}  // namespace hictk
