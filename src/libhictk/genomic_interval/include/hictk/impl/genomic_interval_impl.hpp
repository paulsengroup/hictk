// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string>
#include <string_view>

#include "hictk/chromosome.hpp"
#include "hictk/hash.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/reference.hpp"

namespace hictk {

inline GenomicInterval::GenomicInterval(const Chromosome &chrom_) noexcept
    : GenomicInterval(chrom_, 0, chrom_.size()) {}

inline GenomicInterval::GenomicInterval(const Chromosome &chrom_, std::uint32_t start_,
                                        std::uint32_t end_) noexcept
    : _chrom(chrom_), _start(start_), _end(end_) {
  assert(_start <= _end);
}

inline GenomicInterval::operator bool() const noexcept { return !!chrom(); }

inline bool GenomicInterval::operator==(const GenomicInterval &other) const noexcept {
  // clang-format off
  return chrom() == other.chrom() &&
         start() == other.start() &&
         end() == other.end();
  // clang-format on
}
inline bool GenomicInterval::operator!=(const GenomicInterval &other) const noexcept {
  return !(*this == other);
}

inline bool GenomicInterval::operator<(const GenomicInterval &other) const noexcept {
  if (chrom() != other.chrom()) {
    return chrom() < other.chrom();
  }

  if (start() != other.start()) {
    return start() < other.start();
  }

  return end() < other.end();
}

inline bool GenomicInterval::operator<=(const GenomicInterval &other) const noexcept {
  if (chrom() != other.chrom()) {
    return chrom() <= other.chrom();
  }

  if (start() != other.start()) {
    return start() <= other.start();
  }

  return end() <= other.end();
}

inline bool GenomicInterval::operator>(const GenomicInterval &other) const noexcept {
  if (chrom() != other.chrom()) {
    return chrom() > other.chrom();
  }

  if (start() != other.start()) {
    return start() > other.start();
  }

  return end() > other.end();
}

inline bool GenomicInterval::operator>=(const GenomicInterval &other) const noexcept {
  if (chrom() != other.chrom()) {
    return chrom() >= other.chrom();
  }

  if (start() != other.start()) {
    return start() >= other.start();
  }

  return end() >= other.end();
}

inline const Chromosome &GenomicInterval::chrom() const noexcept { return _chrom; }
constexpr std::uint32_t GenomicInterval::start() const noexcept { return _start; }
constexpr std::uint32_t GenomicInterval::end() const noexcept { return _end; }
constexpr std::uint32_t GenomicInterval::size() const noexcept { return _end - _start; }

inline GenomicInterval GenomicInterval::parse(const Reference &chroms, std::string query,
                                              Type type) {
  if (type == Type::UCSC) {
    return GenomicInterval::parse_ucsc(chroms, query);
  }
  return GenomicInterval::parse_bed(chroms, query);
}

inline GenomicInterval GenomicInterval::parse_ucsc(const Reference &chroms, std::string buffer) {
  if (buffer.empty()) {
    throw std::runtime_error("query is empty");
  }

  if (buffer.back() == '\r') {
    buffer.resize(buffer.size() - 1);
  }

  if (const auto match = chroms.find(buffer); match != chroms.end()) {
    return GenomicInterval{*match};
  }

  const auto p1 = buffer.find_last_of(':');
  auto p2 = buffer.find_last_of('-');

  if (p1 == std::string::npos && p2 == std::string::npos) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid chromosome \"{0}\" in query \"{0}\""), buffer));
  }

  if (p1 == std::string::npos || p2 == std::string::npos || p1 > p2) {
    throw std::runtime_error(fmt::format(FMT_STRING("query \"{}\" is malformed"), buffer));
  }

  if (buffer.find(',', p1) != std::string::npos) {
    buffer.erase(std::remove(buffer.begin() + std::ptrdiff_t(p1), buffer.end(), ','), buffer.end());
    p2 = buffer.find_last_of('-');
  }

  buffer[p1] = '\t';
  buffer[p2] = '\t';

  return GenomicInterval::parse_bed(chroms, buffer);
}

inline GenomicInterval GenomicInterval::parse_bed(const Reference &chroms, std::string_view buffer,
                                                  char sep) {
  if (buffer.empty()) {
    throw std::runtime_error("interval is empty");
  }

  if (buffer.back() == '\r') {
    buffer = buffer.substr(0, buffer.size() - 1);
  }

  const auto p1 = buffer.find(sep);
  const auto p2 = buffer.find(sep, p1 + 1);

  if (p1 == std::string_view::npos || p2 == std::string_view::npos || p1 > p2) {
    throw std::runtime_error(fmt::format(FMT_STRING("interval \"{}\" is malformed"), buffer));
  }

  const auto chrom_name = buffer.substr(0, p1);
  const auto start_pos_str = buffer.substr(p1 + 1, p2 - (p1 + 1));
  const auto end_pos_str = buffer.substr(p2 + 1);

  const auto match = chroms.find(chrom_name);
  if (match == chroms.end()) {
    throw std::runtime_error(fmt::format(FMT_STRING("invalid chromosome \"{}\" in interval \"{}\""),
                                         chrom_name, buffer));
  }

  if (start_pos_str.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("interval \"{}\" is malformed: missing start position"), buffer));
  }

  if (end_pos_str.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("interval \"{}\" is malformed: missing end position"), buffer));
  }

  GenomicInterval gi{*match};

  try {
    internal::parse_numeric_or_throw(start_pos_str, gi._start);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid start position \"{}\" in interval \"{}\": {}"),
                    start_pos_str, buffer, e.what()));
  }
  try {
    internal::parse_numeric_or_throw(end_pos_str, gi._end);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid end position \"{}\" in interval \"{}\": {}"), end_pos_str,
                    buffer, e.what()));
  }

  if (gi._end > gi.chrom().size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid end position \"{0}\" in interval \"{1}\": end position is "
                               "greater than the chromosome size ({0} > {2})"),
                    gi._end, buffer, gi.chrom().size()));
  }

  if (gi._start >= gi._end) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid interval \"{}\": query end position should be "
                               "greater than the start position ({} >= {})"),
                    buffer, gi._start, gi._end));
  }

  return gi;
}

}  // namespace hictk

inline std::size_t std::hash<hictk::GenomicInterval>::operator()(
    const hictk::GenomicInterval &gi) const {
  return hictk::internal::hash_combine(0, gi.chrom(), gi.start(), gi.end());
}
