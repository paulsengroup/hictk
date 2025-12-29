// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/bin.hpp"

#include <cstdint>
#include <utility>

#include "hictk/chromosome.hpp"
#include "hictk/genomic_interval.hpp"

namespace hictk {

Bin::Bin(Chromosome chrom_, std::uint32_t start_, std::uint32_t end_) noexcept
    : Bin(null_id, rel_null_id, std::move(chrom_), start_, end_) {}

Bin::Bin(std::uint64_t id_, std::uint32_t rel_id_, Chromosome chrom_, std::uint32_t start_,
         std::uint32_t end_) noexcept
    : Bin(id_, rel_id_, {std::move(chrom_), start_, end_}) {}

Bin::Bin(GenomicInterval interval) noexcept : Bin(null_id, rel_null_id, std::move(interval)) {}

Bin::Bin(std::uint64_t id_, std::uint32_t rel_id_, GenomicInterval interval) noexcept
    : _id(id_), _rel_id(rel_id_), _interval(std::move(interval)) {}

Bin::operator bool() const noexcept { return !!chrom(); }

bool Bin::operator==(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() == other.id();
  }
  return _interval == other._interval;
}
bool Bin::operator!=(const Bin &other) const noexcept { return !(*this == other); }

bool Bin::operator<(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() < other.id();
  }
  return _interval < other._interval;
}

bool Bin::operator<=(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() <= other.id();
  }
  return _interval <= other._interval;
}

bool Bin::operator>(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() > other.id();
  }
  return _interval > other._interval;
}

bool Bin::operator>=(const Bin &other) const noexcept {
  if (!has_null_id() && !other.has_null_id()) {
    return id() >= other.id();
  }
  return _interval >= other._interval;
}

const GenomicInterval &Bin::interval() const noexcept { return _interval; }
const Chromosome &Bin::chrom() const noexcept { return interval().chrom(); }

}  // namespace hictk

std::size_t std::hash<hictk::Bin>::operator()(const hictk::Bin &b) const noexcept {
  return hictk::internal::hash_combine(0, b.id(), b.interval());
}
