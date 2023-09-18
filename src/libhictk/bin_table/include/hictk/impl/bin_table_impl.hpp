// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
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
#include "hictk/genomic_interval.hpp"
#include "hictk/suppress_warnings.hpp"

namespace hictk {  // NOLINT

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

inline BinTable::BinTable(Reference chroms, std::uint32_t bin_size, std::size_t bin_offset)
    : _chroms(std::move(chroms)),
      _num_bins_prefix_sum(compute_num_bins_prefix_sum(_chroms, bin_size, bin_offset)),
      _bin_size(bin_size) {
  assert(bin_size != 0);
}

template <typename ChromIt>
inline BinTable::BinTable(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size,
                          std::size_t bin_offset)
    : BinTable(Reference(first_chrom, last_chrom), bin_size, bin_offset) {}

template <typename ChromNameIt, typename ChromSizeIt>
inline BinTable::BinTable(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name,
                          ChromSizeIt first_chrom_size, std::uint32_t bin_size,
                          std::size_t bin_offset)
    : BinTable(Reference(first_chrom_name, last_chrom_name, first_chrom_size), bin_size,
               bin_offset) {}

inline std::size_t BinTable::size() const noexcept {
  if (_num_bins_prefix_sum.empty()) {
    return 0;
  }
  return static_cast<std::size_t>(_num_bins_prefix_sum.back() - _num_bins_prefix_sum.front());
}

inline bool BinTable::empty() const noexcept { return size() == 0; }

inline std::size_t BinTable::num_chromosomes() const { return _chroms.size(); }

constexpr std::uint32_t BinTable::bin_size() const noexcept { return _bin_size; }

constexpr const Reference &BinTable::chromosomes() const noexcept { return _chroms; }

constexpr const std::vector<std::uint64_t> &BinTable::num_bin_prefix_sum() const noexcept {
  return _num_bins_prefix_sum;
}

inline auto BinTable::begin() const -> iterator { return iterator(*this); }
inline auto BinTable::end() const -> iterator { return iterator::make_end_iterator(*this); }
inline auto BinTable::cbegin() const -> iterator { return begin(); }
inline auto BinTable::cend() const -> iterator { return end(); }

constexpr std::uint32_t BinTable::iterator::bin_size() const noexcept {
  assert(_bin_table);
  return _bin_table->bin_size();
}

constexpr std::size_t BinTable::iterator::bin_id() const noexcept {
  return _chrom_bin_id + _rel_bin_id;
}

inline BinTableConcrete BinTable::concretize() const {
  std::vector<Chromosome> chroms(size());
  std::vector<std::uint32_t> starts(size());
  std::vector<std::uint32_t> ends(size());

  std::size_t i = 0;
  for (const auto &bin : *this) {
    chroms[i] = bin.chrom();
    starts[i] = bin.start();
    ends[i++] = bin.end();
  }
  assert(i == chroms.size());

  return BinTableConcrete{chroms, starts, ends};
}

inline bool BinTable::operator==(const BinTable &other) const {
  return _bin_size == other._bin_size && _chroms == other._chroms;
}
inline bool BinTable::operator!=(const BinTable &other) const { return !(*this == other); }

inline BinTable BinTable::subset(const Chromosome &chrom) const {
  // GCC8 fails to compile when using if constexpr instead #ifndef
  // See: https://github.com/fmtlib/fmt/issues/1455
#ifndef NDEBUG
  if (!_chroms.contains(chrom)) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom.name()));
  }
#endif
  const auto offset = at(chrom, 0).id();
  return {Reference{chrom}, _bin_size, offset};
}
inline BinTable BinTable::subset(std::string_view chrom_name) const {
  return subset(_chroms.at(chrom_name));
}
inline BinTable BinTable::subset(std::uint32_t chrom_id) const {
  return subset(_chroms.at(chrom_id));
}

inline auto BinTable::find_overlap(const GenomicInterval &query) const
    -> std::pair<BinTable::iterator, BinTable::iterator> {
  return find_overlap(query.chrom(), query.start(), query.end());
}

inline auto BinTable::find_overlap(const Chromosome &chrom, std::uint32_t start,
                                   std::uint32_t end) const
    -> std::pair<BinTable::iterator, BinTable::iterator> {
  assert(start < end);

  const auto bin1_id = at(chrom, start).id();
  const auto bin2_id = at(chrom, end - (std::min)(end, 1U)).id();

  return std::make_pair(begin() + bin1_id, begin() + bin2_id + 1);
}
inline auto BinTable::find_overlap(std::string_view chrom_name, std::uint32_t start,
                                   std::uint32_t end) const
    -> std::pair<BinTable::iterator, BinTable::iterator> {
  return find_overlap(_chroms.at(chrom_name), start, end);
}
inline auto BinTable::find_overlap(std::uint32_t chrom_id, std::uint32_t start,
                                   std::uint32_t end) const
    -> std::pair<BinTable::iterator, BinTable::iterator> {
  return find_overlap(_chroms.at(chrom_id), start, end);
}

inline Bin BinTable::at(std::uint64_t bin_id) const {
  // I tried benchmarking linear search as well as std::set (including third-party implementations).
  // Binary search and find on flat vectors are always faster for a reasonable number of chromosomes
  // (e.g. 5-100) and have fairly similar performanc.
  // Linear search is however better in practice because chromosomes are usually sorted by (approx.)
  // size, with unplaced scaffolds etc. ending up last.
  auto match = std::find_if(_num_bins_prefix_sum.begin(), _num_bins_prefix_sum.end(),
                            [&](const auto n) { return n > bin_id; });

  if (match == _num_bins_prefix_sum.end()) {
    throw std::out_of_range(fmt::format(FMT_STRING("bin id {} not found: out of range"), bin_id));
  }
  assert(match != _num_bins_prefix_sum.begin());

  const auto chrom_id =
      static_cast<std::uint32_t>(std::distance(_num_bins_prefix_sum.begin(), --match));
  return at_hint(bin_id, _chroms[chrom_id]);
}

inline Bin BinTable::at_hint(std::uint64_t bin_id, const Chromosome &chrom) const {
  const auto offset = _num_bins_prefix_sum[chrom.id()];
  const auto relative_bin_id = bin_id - offset;
  const auto start = static_cast<uint32_t>(relative_bin_id * bin_size());
  assert(start < chrom.size());
  const auto end = (std::min)(start + bin_size(), chrom.size());

  return {bin_id, static_cast<std::uint32_t>(relative_bin_id), chrom, start, end};
}

inline std::pair<Bin, Bin> BinTable::at(const GenomicInterval &gi) const {
  const auto [bin1_id, bin2_id] = map_to_bin_ids(gi);
  return std::make_pair(at_hint(bin1_id, gi.chrom()), at_hint(bin2_id, gi.chrom()));
}
inline Bin BinTable::at(const Chromosome &chrom, std::uint32_t pos) const {
  return at_hint(map_to_bin_id(chrom, pos), chrom);
}
inline Bin BinTable::at(std::string_view chrom_name, std::uint32_t pos) const {
  return at(map_to_bin_id(chrom_name, pos));
}
inline Bin BinTable::at(std::uint32_t chrom_id, std::uint32_t pos) const {
  return at(map_to_bin_id(chrom_id, pos));
}

inline std::pair<std::uint64_t, std::uint64_t> BinTable::map_to_bin_ids(
    const GenomicInterval &gi) const {
  return std::make_pair(map_to_bin_id(gi.chrom(), gi.start()),
                        map_to_bin_id(gi.chrom(), gi.end() - (std::min)(gi.end(), 1U)));
}

inline std::uint64_t BinTable::map_to_bin_id(const Chromosome &chrom, std::uint32_t pos) const {
  // GCC8 fails to compile when using if constexpr instead #ifndef
  // See: https://github.com/fmtlib/fmt/issues/1455
#ifndef NDEBUG
  if (!_chroms.contains(chrom)) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom.name()));
  }
#endif

  if (pos > chrom.size()) {
    throw std::out_of_range(fmt::format(
        FMT_STRING("position is greater than chromosome size: {} > {}"), pos, chrom.size()));
  }

  const auto bin_offset = _num_bins_prefix_sum[chrom.id()] - _num_bins_prefix_sum.front();

  return bin_offset + static_cast<std::uint64_t>(pos / bin_size());
}

inline std::uint64_t BinTable::map_to_bin_id(std::string_view chrom_name, std::uint32_t pos) const {
  return map_to_bin_id(_chroms.at(chrom_name), pos);
}

inline std::uint64_t BinTable::map_to_bin_id(std::uint32_t chrom_id, std::uint32_t pos) const {
  return map_to_bin_id(_chroms.at(chrom_id), pos);
}

inline std::vector<std::uint64_t> BinTable::compute_num_bins_prefix_sum(const Reference &chroms,
                                                                        std::uint32_t bin_size,
                                                                        std::size_t bin_offset) {
  assert(bin_size != 0);

  DISABLE_WARNING_PUSH
  DISABLE_WARNING_NULL_DEREF
  std::vector<std::uint64_t> prefix_sum(chroms.size() + 1);
  prefix_sum.front() = bin_offset;
  DISABLE_WARNING_POP

  // I am using transform instead of inclusive_scan because the latter is not always available
  std::transform(chroms.begin(), chroms.end(), prefix_sum.begin() + 1,
                 [&, sum = bin_offset](const Chromosome &chrom) mutable {
                   if (chrom.is_all()) {
                     return sum;
                   }
                   const auto num_bins = (chrom.size() + bin_size - 1) / bin_size;
                   return sum += static_cast<std::uint64_t>(num_bins);
                 });

  return prefix_sum;
}

inline BinTable::iterator::iterator(const BinTable &bin_table) noexcept
    : _bin_table{&bin_table}, _chrom_bin_id(_bin_table->_num_bins_prefix_sum.front()) {
  if (_bin_table->chromosomes().at(_chrom_id).is_all()) {
    _chrom_id++;
  }
}

constexpr bool BinTable::iterator::operator==(const iterator &other) const noexcept {
  // clang-format off
  return _bin_table == other._bin_table &&
         _chrom_id == other._chrom_id &&
         _rel_bin_id == other._rel_bin_id;
  // clang-format on
}

constexpr bool BinTable::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

constexpr bool BinTable::iterator::operator<(const iterator &other) const noexcept {
  return bin_id() < other.bin_id();
}

constexpr bool BinTable::iterator::operator<=(const iterator &other) const noexcept {
  return bin_id() <= other.bin_id();
}

constexpr bool BinTable::iterator::operator>(const iterator &other) const noexcept {
  return bin_id() > other.bin_id();
}

constexpr bool BinTable::iterator::operator>=(const iterator &other) const noexcept {
  return bin_id() >= other.bin_id();
}

inline auto BinTable::iterator::make_end_iterator(const BinTable &table) noexcept -> iterator {
  iterator it(table);

  it._chrom_id = nchrom;
  it._rel_bin_id = null_rel_bin_id;
  return it;
}

inline auto BinTable::iterator::operator*() const -> value_type {
  assert(_bin_table);

  const auto &chrom = chromosome();
  const auto bin_size = this->bin_size();

  const auto start = std::min(_rel_bin_id * bin_size, chrom.size());
  const auto end = std::min(start + bin_size, chrom.size());

  return value_type{bin_id(), _rel_bin_id, chrom, start, end};
}

inline auto BinTable::iterator::operator++() -> iterator & {
  assert(_bin_table);
  if (_chrom_id == nchrom) {
    return *this;
  }

  if (++_rel_bin_id >= compute_num_chrom_bins()) {
    if (_chrom_id + 1 >= num_chromosomes()) {
      return *this = make_end_iterator(*_bin_table);
    }
    ++_chrom_id;
    _chrom_bin_id = compute_bin_offset();
    _rel_bin_id = 0;
  }

  return *this;
}

inline auto BinTable::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

inline auto BinTable::iterator::operator+=(std::size_t i) -> iterator & {
  assert(_bin_table);
  if (_chrom_id == nchrom) {
    if (i == 0) {
      return *this;
    }
    throw std::out_of_range("BinTable::iterator: caught attempt to increment iterator past end()");
  }

  const auto ii = static_cast<std::uint32_t>(i);
  const auto num_bins = compute_num_chrom_bins();
  if (_rel_bin_id + ii < num_bins) {
    _rel_bin_id += ii;
    return *this;
  }

  _chrom_id++;
  _chrom_bin_id = compute_bin_offset();
  i -= (num_bins - _rel_bin_id);
  _rel_bin_id = 0;
  return *this += i;
}

inline auto BinTable::iterator::operator+(std::size_t i) const -> iterator {
  auto it = *this;
  return it += i;
}

inline auto BinTable::iterator::operator--() -> iterator & {
  assert(_bin_table);
  if (bin_id() == 0) {
    return *this;
  }

  if (_rel_bin_id == null_rel_bin_id) {
    assert(*this == make_end_iterator(*_bin_table));
    _chrom_id = static_cast<std::uint32_t>(num_chromosomes() - 1);
    _chrom_bin_id = compute_bin_offset();
    _rel_bin_id = compute_num_chrom_bins() - 1;
    return *this;
  }

  if (_rel_bin_id-- == 0) {
    _chrom_id--;
    _chrom_bin_id = compute_bin_offset();
    _rel_bin_id = compute_num_chrom_bins() - 1;
  }

  return *this;
}

inline auto BinTable::iterator::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  return it;
}

inline auto BinTable::iterator::operator-=(std::size_t i) -> iterator & {
  assert(_bin_table);

  if (_chrom_id == 0 && _rel_bin_id == 0 && i == 0) {
    return *this;
  }

  if (_chrom_id == 0 && _rel_bin_id < i) {
    throw std::out_of_range(
        "BinTable::iterator: caught attempt to decrement iterator past begin()");
  }

  if (_rel_bin_id == null_rel_bin_id) {
    assert(*this == make_end_iterator(*_bin_table));
    _chrom_id = static_cast<std::uint32_t>(num_chromosomes() - 1);
    _chrom_bin_id = compute_bin_offset();
    _rel_bin_id = compute_num_chrom_bins();
    return *this -= i;
  }

  if (i <= _rel_bin_id) {
    _rel_bin_id -= static_cast<std::uint32_t>(i);
    return *this;
  }

  _chrom_id--;
  _chrom_bin_id = compute_bin_offset();
  i -= _rel_bin_id;
  _rel_bin_id = compute_num_chrom_bins();
  return *this -= i;
}

inline auto BinTable::iterator::operator-(std::size_t i) const -> iterator {
  auto it = *this;
  return it -= i;
}

inline auto BinTable::iterator::operator-(const iterator &other) const -> difference_type {
  assert(_bin_table);
  assert(other._bin_table);

  const auto offset1 = _chrom_id == nchrom ? _bin_table->size() : bin_id();
  const auto offset2 = other._chrom_id == nchrom ? other._bin_table->size() : other.bin_id();

  return static_cast<difference_type>(offset1) - static_cast<difference_type>(offset2);
}

inline auto BinTable::iterator::operator[](std::size_t i) const -> iterator { return (*this + i); }

inline const Chromosome &BinTable::iterator::chromosome() const { return chromosome(_chrom_id); }

inline const Chromosome &BinTable::iterator::chromosome(std::uint32_t chrom_id) const {
  return _bin_table->chromosomes().at(chrom_id);
}

inline std::uint32_t BinTable::iterator::compute_num_chrom_bins() const noexcept {
  assert(_bin_table);

  const auto chrom_size = chromosome().size();
  const auto bin_size = this->bin_size();

  return (chrom_size + bin_size - 1) / bin_size;
}

inline std::size_t BinTable::iterator::compute_bin_offset() const noexcept {
  return _bin_table->at(_chrom_id, 0).id();
}

inline std::size_t BinTable::iterator::num_chromosomes() const noexcept {
  assert(_bin_table);

  return _bin_table->num_chromosomes();
}

}  // namespace hictk

inline std::size_t std::hash<hictk::Bin>::operator()(const hictk::Bin &b) const {
  return hictk::internal::hash_combine(0, b.id(), b.interval());
}
