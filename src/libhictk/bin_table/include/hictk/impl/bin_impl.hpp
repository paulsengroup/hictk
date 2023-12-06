// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>

#include "hictk/chromosome.hpp"
#include "hictk/genomic_interval.hpp"

namespace hictk {

inline Bin::Bin(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end_) noexcept
    : Bin(Bin::null_id, Bin::rel_null_id, chrom_, start_, end_) {}

inline Bin::Bin(std::uint64_t id_, std::uint32_t rel_id_, const Chromosome &chrom_,
                std::uint32_t start_, std::uint32_t end_) noexcept
    : _id(id_), _rel_id(rel_id_), _interval(chrom_, start_, end_) {}

inline Bin::Bin(GenomicInterval interval) noexcept
    : Bin(Bin::null_id, Bin::rel_null_id, std::move(interval)) {}

inline Bin::Bin(std::uint64_t id_, std::uint32_t rel_id_, GenomicInterval interval) noexcept
    : _id(id_), _rel_id(rel_id_), _interval(std::move(interval)) {}

inline Bin::operator bool() const noexcept { return !!chrom(); }

inline bool Bin::operator==(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() == other.id();
  }
  return _interval == other._interval;
}
inline bool Bin::operator!=(const Bin &other) const noexcept { return !(*this == other); }

inline bool Bin::operator<(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() < other.id();
  }
  return _interval < other._interval;
}

inline bool Bin::operator<=(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() <= other.id();
  }
  return _interval <= other._interval;
}

inline bool Bin::operator>(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() > other.id();
  }
  return _interval > other._interval;
}

inline bool Bin::operator>=(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() >= other.id();
  }
  return _interval >= other._interval;
}

constexpr std::uint64_t Bin::id() const noexcept { return _id; }
constexpr std::uint32_t Bin::rel_id() const noexcept { return _rel_id; }
inline const GenomicInterval &Bin::interval() const noexcept { return _interval; }
inline const Chromosome &Bin::chrom() const noexcept { return interval().chrom(); }
constexpr std::uint32_t Bin::start() const noexcept { return _interval.start(); }
constexpr std::uint32_t Bin::end() const noexcept { return _interval.end(); }

constexpr bool Bin::has_null_id() const noexcept { return id() == Bin::null_id; }

};  // namespace hictk

inline std::size_t std::hash<hictk::Bin>::operator()(const hictk::Bin &b) const {
  return hictk::internal::hash_combine(0, b.id(), b.interval());
}
