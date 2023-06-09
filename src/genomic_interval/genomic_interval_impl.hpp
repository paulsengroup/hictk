// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/numeric_utils.hpp"

namespace hictk {

inline GenomicInterval::GenomicInterval(const Chromosome &chrom_) noexcept
    : GenomicInterval(chrom_, 0, chrom_.size()) {}

inline GenomicInterval::GenomicInterval(const Chromosome &chrom_, std::uint32_t start_,
                                        std::uint32_t end_) noexcept
    : _chrom(chrom_), _start(start_), _end(end_) {
  assert(_start <= _end);
}

inline GenomicInterval::operator bool() const noexcept { return !!this->chrom(); }

inline bool GenomicInterval::operator==(const GenomicInterval &other) const noexcept {
  // clang-format off
  return this->chrom() == other.chrom() &&
         this->start() == other.start() &&
         this->end() == other.end();
  // clang-format on
}
inline bool GenomicInterval::operator!=(const GenomicInterval &other) const noexcept {
  return !(*this == other);
}

inline bool GenomicInterval::operator<(const GenomicInterval &other) const noexcept {
  if (this->chrom() != other.chrom()) {
    return this->chrom() < other.chrom();
  }

  if (this->start() != other.start()) {
    return this->start() < other.start();
  }

  return this->end() < other.end();
}

inline bool GenomicInterval::operator<=(const GenomicInterval &other) const noexcept {
  if (this->chrom() != other.chrom()) {
    return this->chrom() <= other.chrom();
  }

  if (this->start() != other.start()) {
    return this->start() <= other.start();
  }

  return this->end() <= other.end();
}

inline bool GenomicInterval::operator>(const GenomicInterval &other) const noexcept {
  if (this->chrom() != other.chrom()) {
    return this->chrom() > other.chrom();
  }

  if (this->start() != other.start()) {
    return this->start() > other.start();
  }

  return this->end() > other.end();
}

inline bool GenomicInterval::operator>=(const GenomicInterval &other) const noexcept {
  if (this->chrom() != other.chrom()) {
    return this->chrom() >= other.chrom();
  }

  if (this->start() != other.start()) {
    return this->start() >= other.start();
  }

  return this->end() >= other.end();
}

inline const Chromosome &GenomicInterval::chrom() const noexcept { return this->_chrom; }
constexpr std::uint32_t GenomicInterval::start() const noexcept { return this->_start; }
constexpr std::uint32_t GenomicInterval::end() const noexcept { return this->_end; }

inline GenomicInterval GenomicInterval::parse(const Reference &chroms, std::string query,
                                              Type type) {
  if (type == Type::UCSC) {
    return GenomicInterval::parse_ucsc(chroms, query);
  }
  return GenomicInterval::parse_bed(chroms, query);
}

inline GenomicInterval GenomicInterval::parse_ucsc(const Reference &chroms, std::string query) {
  if (query.empty()) {
    throw std::runtime_error("query is empty");
  }

  if (const auto match = chroms.find(query); match != chroms.end()) {
    return GenomicInterval{*match};
  }

  const auto p1 = query.find_last_of(':');
  const auto p2 = query.find_last_of('-');

  if (p1 == std::string::npos && p2 == std::string::npos) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid chromosome \"{0}\" in query \"{0}\""), query));
  }

  if (p1 == std::string::npos || p2 == std::string::npos || p1 > p2) {
    throw std::runtime_error(fmt::format(FMT_STRING("query \"{}\" is malformed"), query));
  }

  if (query.find(',', p1) != std::string::npos) {
    query.erase(std::remove(query.begin() + std::ptrdiff_t(p1), query.end(), ','));
  }

  query[p1] = '\t';
  query[p2] = '\t';

  return GenomicInterval::parse_bed(chroms, query);
}

inline GenomicInterval GenomicInterval::parse_bed(const Reference &chroms, std::string_view query,
                                                  char sep) {
  if (query.empty()) {
    throw std::runtime_error("query is empty");
  }

  const auto p1 = query.find(sep);
  const auto p2 = query.find(sep, p1 + 1);

  if (p1 == std::string_view::npos || p2 == std::string_view::npos || p1 > p2) {
    throw std::runtime_error(fmt::format(FMT_STRING("query \"{}\" is malformed"), query));
  }

  const auto chrom_name = query.substr(0, p1);
  const auto start_pos_str = query.substr(p1 + 1, p2 - (p1 + 1));
  const auto end_pos_str = query.substr(p2 + 1);

  const auto match = chroms.find(chrom_name);
  if (match == chroms.end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid chromosome \"{}\" in query \"{}\""), chrom_name, query));
  }

  if (start_pos_str.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("query \"{}\" is malformed: missing start position"), query));
  }

  if (end_pos_str.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("query \"{}\" is malformed: missing end position"), query));
  }

  GenomicInterval gi{*match};

  try {
    internal::parse_numeric_or_throw(start_pos_str, gi._start);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid start position \"{}\" in query \"{}\": {}"), start_pos_str,
                    query, e.what()));
  }
  try {
    internal::parse_numeric_or_throw(end_pos_str, gi._end);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid end position \"{}\" in query \"{}\": {}"), end_pos_str,
                    query, e.what()));
  }

  if (gi._end > gi.chrom().size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid end position \"{0}\" in query \"{1}\": end position is "
                               "greater than the chromosome size ({0} > {2})"),
                    gi._end, query, gi.chrom().size()));
  }

  if (gi._start >= gi._end) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("invalid query \"{}\": query end position should be "
                               "greater than the start position ({} >= {})"),
                    query, gi._start, gi._end));
  }

  return gi;
}

}  // namespace hictk

inline std::size_t std::hash<hictk::GenomicInterval>::operator()(
    const hictk::GenomicInterval &gi) const {
  return hictk::internal::hash_combine(0, gi.chrom(), gi.start(), gi.end());
}
