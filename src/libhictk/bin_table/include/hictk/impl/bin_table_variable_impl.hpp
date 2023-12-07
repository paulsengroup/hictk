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

template <typename I>
inline BinTableVariable<I>::BinTableVariable(Reference chroms, const std::vector<I> &start_pos,
                                             const std::vector<I> &end_pos, I bin_offset)
    : _chroms(std::move(chroms)) {
  assert(start_pos.size() == end_pos.size());
  _bin_end_prefix_sum[0] = bin_offset;
  if (start_pos.empty()) {
    return;
  }

  _bin_end_prefix_sum.reserve(start_pos.size() + 1);
  _bin_end_prefix_sum.push_back(end_pos.front());

  for (std::size_t i = 1; i < start_pos.size(); ++i) {
    if (start_pos[i] <= start_pos[i - 1]) {
      _num_bins_prefix_sum.push_back(_bin_end_prefix_sum.size() - 1);
    }
    _bin_end_prefix_sum.push_back(_bin_end_prefix_sum.back() + end_pos[i] - start_pos[i]);
  }
  _num_bins_prefix_sum.push_back(_bin_end_prefix_sum.size() - 1);
}

template <typename I>
inline std::size_t BinTableVariable<I>::size() const noexcept {
  return _bin_end_prefix_sum.size() - 1;
}

template <typename I>
inline bool BinTableVariable<I>::empty() const noexcept {
  return size() == 0;
}

template <typename I>
inline std::size_t BinTableVariable<I>::num_chromosomes() const {
  return chromosomes().size();
}

template <typename I>
constexpr const Reference &BinTableVariable<I>::chromosomes() const noexcept {
  return _chroms;
}

template <typename I>
constexpr const std::vector<std::uint64_t> &BinTableVariable<I>::num_bin_prefix_sum()
    const noexcept {
  return _num_bins_prefix_sum;
}

template <typename I>
inline auto BinTableVariable<I>::begin() const -> iterator {
  volatile auto chrom = *chromosomes().begin();
  return iterator(*this);
}

template <typename I>
inline auto BinTableVariable<I>::end() const -> iterator {
  return iterator::make_end_iterator(*this);
}

template <typename I>
inline auto BinTableVariable<I>::cbegin() const -> iterator {
  return begin();
}

template <typename I>
inline auto BinTableVariable<I>::cend() const -> iterator {
  return end();
}

template <typename I>
inline bool BinTableVariable<I>::operator==(const BinTableVariable &other) const {
  return _chroms == other._chroms &&
         std::equal(_num_bins_prefix_sum.begin(), _num_bins_prefix_sum.end(),
                    other._num_bins_prefix_sum.begin());
}

template <typename I>
inline bool BinTableVariable<I>::operator!=(const BinTableVariable &other) const {
  return !(*this == other);
}

template <typename I>
inline BinTableVariable<I> BinTableVariable<I>::subset(const Chromosome &chrom) const {
  // GCC8 fails to compile when using if constexpr instead #ifndef
  // See: https://github.com/fmtlib/fmt/issues/1455
#ifndef NDEBUG
  if (!_chroms.contains(chrom)) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom.name()));
  }
#endif
  std::vector<I> start_pos{};
  std::vector<I> end_pos{};

  const auto [first_bin, last_bin] = find_overlap(chrom, 0, chrom.size());
  std::for_each(first_bin, last_bin, [&](const Bin &b) {
    start_pos.push_back(b.start());
    end_pos.push_back(b.end());
  });

  return {Reference{chrom}, start_pos, end_pos, start_pos.front()};
}

template <typename I>
inline BinTableVariable<I> BinTableVariable<I>::subset(std::string_view chrom_name) const {
  return subset(_chroms.at(chrom_name));
}

template <typename I>
inline BinTableVariable<I> BinTableVariable<I>::subset(std::uint32_t chrom_id) const {
  return subset(_chroms.at(chrom_id));
}

template <typename I>
inline auto BinTableVariable<I>::find_overlap(const GenomicInterval &query) const
    -> std::pair<BinTableVariable<I>::iterator, BinTableVariable<I>::iterator> {
  return find_overlap(query.chrom(), query.start(), query.end());
}

template <typename I>
inline auto BinTableVariable<I>::find_overlap(const Chromosome &chrom, std::uint32_t start,
                                              std::uint32_t end) const
    -> std::pair<BinTableVariable<I>::iterator, BinTableVariable<I>::iterator> {
  assert(start < end);

  const auto bin1_id = at(chrom, start).id();
  const auto bin2_id = at(chrom, end - (std::min)(end, 1U)).id();

  return std::make_pair(begin() + bin1_id, begin() + bin2_id + 1);
}

template <typename I>
inline auto BinTableVariable<I>::find_overlap(std::string_view chrom_name, std::uint32_t start,
                                              std::uint32_t end) const
    -> std::pair<BinTableVariable<I>::iterator, BinTableVariable<I>::iterator> {
  return find_overlap(_chroms.at(chrom_name), start, end);
}

template <typename I>
inline auto BinTableVariable<I>::find_overlap(std::uint32_t chrom_id, std::uint32_t start,
                                              std::uint32_t end) const
    -> std::pair<BinTableVariable<I>::iterator, BinTableVariable<I>::iterator> {
  return find_overlap(_chroms.at(chrom_id), start, end);
}

template <typename I>
Bin BinTableVariable<I>::at(std::uint64_t bin_id) const {
  // I tried benchmarking linear search as well as std::set (including third-party implementations).
  // Binary search and find on flat vectors are always faster for a reasonable number of chromosomes
  // (e.g. 5-100) and have fairly similar performance.
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

template <typename I>
inline Bin BinTableVariable<I>::at_hint(std::uint64_t bin_id, const Chromosome &chrom) const {
  if (_bin_end_prefix_sum.size() <= bin_id) {
    throw std::out_of_range(fmt::format(FMT_STRING("bin id {} not found: out of range"), bin_id));
  }

  const auto bin_id_offset = _num_bins_prefix_sum[chrom.id()];
  const auto chrom_size_offset = _chroms.chrom_size_prefix_sum()[chrom.id()];
  const auto relative_bin_id = bin_id - bin_id_offset;

  const auto start =
      conditional_static_cast<std::uint32_t>(_bin_end_prefix_sum[bin_id] - chrom_size_offset);
  const auto end =
      conditional_static_cast<std::uint32_t>(_bin_end_prefix_sum[bin_id + 1] - chrom_size_offset);

  if (_bin_end_prefix_sum[bin_id] < chrom_size_offset || end > chrom.size()) {
    throw std::out_of_range(fmt::format(
        FMT_STRING("bin id {} not found using {} as hint: out of range"), bin_id, chrom.name()));
  }

  return {bin_id, static_cast<std::uint32_t>(relative_bin_id), chrom, start, end};
}

template <typename I>
inline std::pair<Bin, Bin> BinTableVariable<I>::at(const GenomicInterval &gi) const {
  const auto [bin1_id, bin2_id] = map_to_bin_ids(gi);
  return std::make_pair(at_hint(bin1_id, gi.chrom()), at_hint(bin2_id, gi.chrom()));
}

template <typename I>
inline Bin BinTableVariable<I>::at(const Chromosome &chrom, std::uint32_t pos) const {
  return at_hint(map_to_bin_id(chrom, pos), chrom);
}

template <typename I>
inline Bin BinTableVariable<I>::at(std::string_view chrom_name, std::uint32_t pos) const {
  return at(map_to_bin_id(chrom_name, pos));
}

template <typename I>
inline Bin BinTableVariable<I>::at(std::uint32_t chrom_id, std::uint32_t pos) const {
  return at(map_to_bin_id(chrom_id, pos));
}

template <typename I>
inline std::pair<std::uint64_t, std::uint64_t> BinTableVariable<I>::map_to_bin_ids(
    const GenomicInterval &gi) const {
  return std::make_pair(map_to_bin_id(gi.chrom(), gi.start()),
                        map_to_bin_id(gi.chrom(), gi.end() - (std::min)(gi.end(), 1U)));
}

template <typename I>
inline std::uint64_t BinTableVariable<I>::map_to_bin_id(const Chromosome &chrom,
                                                        std::uint32_t pos) const {
  // GCC8 fails to compile when using if constexpr instead #ifndef
  // See: https://github.com/fmtlib/fmt/issues/1455
#ifndef NDEBUG
  if (!_chroms.contains(chrom)) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom.name()));
  }
#endif

  if (pos >= chrom.size()) {
    throw std::out_of_range(fmt::format(
        FMT_STRING("position is greater than chromosome size: {} >= {}"), pos, chrom.size()));
  }

  const auto pos_offset = _chroms.chrom_size_prefix_sum()[chrom.id()];
  auto match =
      std::upper_bound(_bin_end_prefix_sum.begin(), _bin_end_prefix_sum.end(), pos + pos_offset);

  return static_cast<std::uint64_t>(std::distance(_bin_end_prefix_sum.begin(), --match));
}

template <typename I>
inline std::uint64_t BinTableVariable<I>::map_to_bin_id(std::string_view chrom_name,
                                                        std::uint32_t pos) const {
  return map_to_bin_id(_chroms.at(chrom_name), pos);
}

template <typename I>
inline std::uint64_t BinTableVariable<I>::map_to_bin_id(std::uint32_t chrom_id,
                                                        std::uint32_t pos) const {
  return map_to_bin_id(_chroms.at(chrom_id), pos);
}

template <typename I>
inline BinTableVariable<I>::iterator::iterator(const BinTableVariable &bin_table) noexcept
    : _bin_table{&bin_table} {
  if (_bin_table->chromosomes().at(_chrom_id).is_all()) {
    _chrom_id++;
  }
  _value = get_bin();
}

template <typename I>
constexpr bool BinTableVariable<I>::iterator::operator==(const iterator &other) const noexcept {
  // clang-format off
  return _bin_table == other._bin_table &&
         _chrom_id == other._chrom_id &&
         _bin_id == other._bin_id;
  // clang-format on
}

template <typename I>
constexpr bool BinTableVariable<I>::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename I>
constexpr bool BinTableVariable<I>::iterator::operator<(const iterator &other) const noexcept {
  return _bin_id < other._bin_id;
}

template <typename I>
constexpr bool BinTableVariable<I>::iterator::operator<=(const iterator &other) const noexcept {
  return _bin_id <= other._bin_id;
}

template <typename I>
constexpr bool BinTableVariable<I>::iterator::operator>(const iterator &other) const noexcept {
  return _bin_id > other._bin_id;
}

template <typename I>
constexpr bool BinTableVariable<I>::iterator::operator>=(const iterator &other) const noexcept {
  return _bin_id >= other._bin_id;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::make_end_iterator(const BinTableVariable &table) noexcept
    -> iterator {
  iterator it(table);

  it._chrom_id = nchrom;
  it._bin_id = table.size();
  return it;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator*() const noexcept -> value_type {
  return _value;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator++() -> iterator & {
  assert(_bin_table);
  if (++_bin_id == _bin_table->size()) {
    *this = make_end_iterator(*_bin_table);
    return *this;
  }
  try {
    _value = get_bin();
    _chrom_id = _value.chrom().id();
  } catch (const std::out_of_range &) {
    throw std::out_of_range(
        "BinTableVariable<I>::iterator: caught attempt to increment iterator past end()");
  }

  return *this;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator+=(std::size_t i) -> iterator & {
  if (i == 0) {
    return *this;
  }
  assert(_bin_table);
  try {
    if (_bin_id > _bin_table->size() - i) {
      throw std::out_of_range("");
    }

    _bin_id += i;
    _value = get_bin();
    _chrom_id = _value.chrom().id();
  } catch (const std::out_of_range &) {
    throw std::out_of_range(
        "BinTableVariable<I>::iterator: caught attempt to increment iterator past end()");
  }
  return *this;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator+(std::size_t i) const -> iterator {
  auto it = *this;
  return it += i;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator--() -> iterator & {
  assert(_bin_table);
  --_bin_id;

  try {
    _value = get_bin();
    _chrom_id = _value.chrom().id();
  } catch (const std::out_of_range &) {
    throw std::out_of_range(
        "BinTableVariable<I>::iterator: caught attempt to decrement iterator past begin()");
  }
  return *this;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator--(int) -> iterator {
  auto it = *this;
  std::ignore = --(*this);
  return it;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator-=(std::size_t i) -> iterator & {
  if (i == 0) {
    return *this;
  }
  assert(_bin_table);
  try {
    if (_bin_id < _bin_table->size() - i) {
      throw std::out_of_range("");
    }
    _bin_id -= i;
    _value = get_bin();
    _chrom_id = _value.chrom().id();
  } catch (const std::out_of_range &) {
    throw std::out_of_range(
        "BinTableVariable<I>::iterator: caught attempt to decrement iterator past begin()");
  }
  return *this;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator-(std::size_t i) const -> iterator {
  auto it = *this;
  return it -= i;
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator-(const iterator &other) const
    -> difference_type {
  assert(_bin_table);
  assert(other._bin_table);

  const auto offset1 = _chrom_id == nchrom ? _bin_table->size() : _bin_id;
  const auto offset2 = other._chrom_id == nchrom ? other._bin_table->size() : other._bin_id;

  return static_cast<difference_type>(offset1) - static_cast<difference_type>(offset2);
}

template <typename I>
inline auto BinTableVariable<I>::iterator::operator[](std::size_t i) const -> iterator {
  return (*this + i);
}

template <typename I>
inline const Chromosome &BinTableVariable<I>::iterator::chromosome() const {
  return chromosome(_chrom_id);
}

template <typename I>
inline const Chromosome &BinTableVariable<I>::iterator::chromosome(std::uint32_t chrom_id) const {
  return _bin_table->chromosomes().at(chrom_id);
}

template <typename I>
inline Bin BinTableVariable<I>::iterator::get_bin() const {
  if (_bin_id == _bin_table->size()) {
    return {};
  }
  try {
    const auto &chrom = chromosome();
    return _bin_table->at_hint(_bin_id, chrom);
  } catch (const std::out_of_range &) {
    return _bin_table->at(_bin_id);
  }
}

}  // namespace hictk
