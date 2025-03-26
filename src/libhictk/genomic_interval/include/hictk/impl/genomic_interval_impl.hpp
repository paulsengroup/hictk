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
#include <limits>
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

inline GenomicInterval::GenomicInterval(Chromosome chrom_, std::uint32_t start_,
                                        std::uint32_t end_) noexcept
    : _chrom(std::move(chrom_)), _start(start_), _end(end_) {
  assert(_start <= _end);
  assert(_end <= _chrom.size());
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
constexpr bool GenomicInterval::empty() const noexcept { return size() == 0; }

inline std::tuple<std::string, std::uint32_t, std::uint32_t> GenomicInterval::parse(
    const std::string &query, Type type) {
  if (type == Type::UCSC) {
    return parse_ucsc(query);
  }
  return parse_bed(query);
}

inline std::tuple<std::string, std::uint32_t, std::uint32_t> GenomicInterval::parse_ucsc(
    std::string buffer) {
  if (buffer.empty()) {
    throw std::runtime_error("query is empty");
  }

  if (buffer.back() == '\r') {
    buffer.pop_back();
  }

  const auto p1 = buffer.find_last_of(':');
  auto p2 = buffer.find_last_of('-');

  if (p1 == std::string::npos && p2 == std::string::npos) {
    return std::make_tuple(std::move(buffer), 0, 0);
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

  return parse_bed(buffer);
}

inline std::tuple<std::string, std::uint32_t, std::uint32_t> GenomicInterval::parse_bed(
    std::string_view buffer, char sep) {
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

  std::string chrom_name{buffer.substr(0, p1)};
  const auto start_pos_str = buffer.substr(p1 + 1, p2 - (p1 + 1));
  const auto end_pos_str = buffer.substr(p2 + 1);

  if (start_pos_str.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("interval \"{}\" is malformed: missing start position"), buffer));
  }

  if (end_pos_str.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("interval \"{}\" is malformed: missing end position"), buffer));
  }

  std::uint32_t start_pos{};
  std::uint32_t end_pos{};

  try {
    internal::parse_numeric_or_throw(start_pos_str, start_pos);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid start position \"{}\" in interval \"{}\": {}"),
                    start_pos_str, buffer, e.what()));
  }
  try {
    internal::parse_numeric_or_throw(end_pos_str, end_pos);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid end position \"{}\" in interval \"{}\": {}"), end_pos_str,
                    buffer, e.what()));
  }

  if (start_pos >= end_pos) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid interval \"{}\": query end position should be "
                               "greater than the start position ({} >= {})"),
                    buffer, start_pos, end_pos));
  }

  return std::make_tuple(std::move(chrom_name), start_pos, end_pos);
}

inline GenomicInterval GenomicInterval::parse(const Reference &chroms, const std::string &query,
                                              Type type) {
  if (type == Type::UCSC) {
    return parse_ucsc(chroms, query);
  }
  return parse_bed(chroms, query);
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

  return parse_bed(chroms, buffer);
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

inline std::uint64_t area(const GenomicInterval &gi, std::uint32_t resolution,
                          bool upper_triangular) noexcept {
  return area(gi, gi, resolution, upper_triangular);
}

inline std::uint64_t area(const GenomicInterval &gi1, const GenomicInterval &gi2,
                          std::uint32_t resolution, bool upper_triangular) noexcept {
  if (gi1.chrom() != gi2.chrom()) {
    return area(gi1.start(), gi1.end(), gi2.start(), gi2.end(), resolution, false);
  }
  return area(gi1.start(), gi1.end(), gi2.start(), gi2.end(), resolution, upper_triangular);
}

inline std::uint64_t area(std::uint32_t start_pos, std::uint32_t end_pos, std::uint32_t resolution,
                          bool upper_triangular) noexcept {
  return area(start_pos, end_pos, start_pos, end_pos, resolution, upper_triangular);
}

inline std::uint64_t area(std::uint32_t start1_pos, std::uint32_t end1_pos,
                          std::uint32_t start2_pos, std::uint32_t end2_pos,
                          std::uint32_t resolution, bool upper_triangular) noexcept {
  assert(start1_pos <= end1_pos);
  assert(start2_pos <= end2_pos);
  assert(resolution > 0);
  if (start1_pos == end1_pos || start2_pos == end2_pos) {
    return 0;
  }

  const auto i1 = static_cast<std::uint64_t>(start1_pos / resolution);
  const auto i2 = static_cast<std::uint64_t>((end1_pos - 1) / resolution);
  const auto j1 = static_cast<std::uint64_t>(start2_pos / resolution);
  const auto j2 = static_cast<std::uint64_t>((end2_pos - 1) / resolution);

  if (i1 == i2 || j1 == j2) {
    return 0;
  }

  if (!upper_triangular) {
    return (i2 - i1) * (j2 - j1);
  }

  if (i2 <= j1 + 1) {
    // area does not cross the matrix diagonal
    return area(start1_pos, end1_pos, start2_pos, end2_pos, resolution, false);
  }

  constexpr auto npos = std::numeric_limits<std::uint64_t>::max();
  const auto intersection_left = j1;
  auto intersection_right = npos;

  if (j1 < i2) {
    for (auto x = j2; x > j1; --x) {
      if (i2 > x) {
        intersection_right = x;
        break;
      }
    }
  }

  // this check is probably redundant
  if (intersection_right == npos) {
    // area does not cross the matrix diagonal
    return area(start1_pos, end1_pos, start2_pos, end2_pos, resolution, false);
  }

  // if we get to this point it means that we are in a scenario like the one depicted below:
  // ------------
  // |          |
  // |    A     |
  // |╲         | <- begin of intersection with the matrix diagonal
  // |  ╲       |
  // |    ╲     |
  // |  B   ╲   |
  // |--------╲ | <- end of intersection with the matrix diagonal
  // |          |
  // |    C     |
  // |          |
  // ------------
  // we want to calculate the area of the trapezoid A.
  // This can be achieved by:
  // - computing the total area delimited by the two genomic intervals
  // - computing the area of triangle B
  // - computing the area of the rectangle C
  // - subtracting B + C from the total area

  const auto total_area = (i2 - i1) * (j2 - j1);
  const auto intersection_width = intersection_right - intersection_left;
  const auto triangle_area = (intersection_width * (intersection_width + 1)) / 2;
  const auto square_area = (intersection_width * (i2 - 1 - intersection_left)) -
                           (intersection_width * intersection_width);

  assert(total_area >= triangle_area + square_area);
  return total_area - triangle_area - square_area;
}

}  // namespace hictk

inline std::size_t std::hash<hictk::GenomicInterval>::operator()(
    const hictk::GenomicInterval &gi) const noexcept {
  return hictk::internal::hash_combine(0, gi.chrom(), gi.start(), gi.end());
}
