// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/bin_table.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <string_view>
#include <utility>
#include <variant>

#include "hictk/bin.hpp"
#include "hictk/bin_table_fixed.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/suppress_warnings.hpp"

namespace hictk {
BinTable::BinTable(Reference chroms, std::uint32_t bin_size, std::size_t bin_offset)
    : BinTable(BinTableFixed(std::move(chroms), bin_size, bin_offset)) {}

// NOLINTNEXTLINE(bugprone-exception-escape)
std::size_t BinTable::size() const noexcept {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.size(); }, _table);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
bool BinTable::empty() const noexcept { return size() == 0; }

// NOLINTNEXTLINE(bugprone-exception-escape)
std::size_t BinTable::num_chromosomes() const noexcept {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.num_chromosomes(); }, _table);
}

auto BinTable::begin() const -> iterator {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return iterator{t.begin()}; }, _table);
}

auto BinTable::end() const -> iterator {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return iterator{t.end()}; }, _table);
}

auto BinTable::cbegin() const -> iterator { return begin(); }

auto BinTable::cend() const -> iterator { return end(); }

BinTable BinTable::subset(const Chromosome &chrom) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return BinTable{t.subset(chrom)}; }, _table);
}

BinTable BinTable::subset(std::string_view chrom_name) const {
  return subset(chromosomes().at(chrom_name));
}

BinTable BinTable::subset(std::uint32_t chrom_id) const {
  return subset(chromosomes().at(chrom_id));
}

auto BinTable::find_overlap(const GenomicInterval &query) const -> std::pair<iterator, iterator> {
  assert(!_table.valueless_by_exception());
  return std::visit(
      [&](const auto &t) {
        auto its = t.find_overlap(query);
        HICTK_DISABLE_WARNING_PUSH
        HICTK_DISABLE_WARNING_MAYBE_UNINITIALIZED
        return std::make_pair(iterator{its.first}, iterator{its.second});
        HICTK_DISABLE_WARNING_POP
      },
      _table);
}

auto BinTable::find_overlap(const Chromosome &chrom, std::uint32_t start, std::uint32_t end) const
    -> std::pair<iterator, iterator> {
  return find_overlap(GenomicInterval{chrom, start, end});
}

auto BinTable::find_overlap(std::string_view chrom_name, std::uint32_t start,
                            std::uint32_t end) const -> std::pair<iterator, iterator> {
  return find_overlap(chromosomes().at(chrom_name), start, end);
}

auto BinTable::find_overlap(std::uint32_t chrom_id, std::uint32_t start, std::uint32_t end) const
    -> std::pair<iterator, iterator> {
  return find_overlap(chromosomes().at(chrom_id), start, end);
}

// Map bin_id to Bin
Bin BinTable::at(std::uint64_t bin_id) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.at(bin_id); }, _table);
}

std::pair<Bin, Bin> BinTable::at(const GenomicInterval &gi) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.at(gi); }, _table);
}

Bin BinTable::at(const Chromosome &chrom, std::uint32_t pos) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.at(chrom, pos); }, _table);
}

Bin BinTable::at(std::string_view chrom_name, std::uint32_t pos) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.at(chrom_name, pos); }, _table);
}

Bin BinTable::at(std::uint32_t chrom_id, std::uint32_t pos) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.at(chrom_id, pos); }, _table);
}

Bin BinTable::at_hint(std::uint64_t bin_id, const Chromosome &chrom) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.at_hint(bin_id, chrom); }, _table);
}

// Map genomic coords to bin_id
std::pair<std::uint64_t, std::uint64_t> BinTable::map_to_bin_ids(const GenomicInterval &gi) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.map_to_bin_ids(gi); }, _table);
}

std::uint64_t BinTable::map_to_bin_id(const Chromosome &chrom, std::uint32_t pos) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.map_to_bin_id(chrom, pos); }, _table);
}

std::uint64_t BinTable::map_to_bin_id(std::string_view chrom_name, std::uint32_t pos) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.map_to_bin_id(chrom_name, pos); }, _table);
}

std::uint64_t BinTable::map_to_bin_id(std::uint32_t chrom_id, std::uint32_t pos) const {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) { return t.map_to_bin_id(chrom_id, pos); }, _table);
}

bool BinTable::operator==(const BinTable &other) const {
  assert(!_table.valueless_by_exception());
  return std::visit(
      [&](const auto &t1) {
        try {
          using BinTableT = remove_cvref_t<decltype(t1)>;
          const auto &t2 = std::get<BinTableT>(other._table);
          return t1 == t2;
        } catch (const std::bad_variant_access &) {
          return false;
        }
      },
      _table);
}

bool BinTable::operator!=(const BinTable &other) const { return !(*this == other); }

auto BinTable::iterator::operator*() const -> value_type {
  assert(!_it.valueless_by_exception());
  return std::visit([&](const auto &it) { return *it; }, _it);
}
auto BinTable::iterator::operator[](difference_type i) const -> value_type {
  assert(!_it.valueless_by_exception());
  return std::visit([&](const auto &it) { return it[i]; }, _it);
}

auto BinTable::iterator::operator++() -> iterator & {
  assert(!_it.valueless_by_exception());
  std::visit([&](auto &it) { std::ignore = ++it; }, _it);
  return *this;
}

auto BinTable::iterator::operator++(int) -> iterator {
  assert(!_it.valueless_by_exception());
  auto old_it = *this;
  std::visit([&](auto &it) { std::ignore = ++it; }, _it);
  return old_it;
}

auto BinTable::iterator::operator+=(difference_type i) -> iterator & {
  assert(!_it.valueless_by_exception());
  std::visit([&](auto &it) { std::ignore = it += i; }, _it);
  return *this;
}

auto BinTable::iterator::operator+(difference_type i) const -> iterator {
  assert(!_it.valueless_by_exception());
  auto it = *this;
  it += i;
  return it;
}

auto BinTable::iterator::operator--() -> iterator & {
  assert(!_it.valueless_by_exception());
  std::visit([&](auto &it) { std::ignore = --it; }, _it);
  return *this;
}

auto BinTable::iterator::operator--(int) -> iterator {
  assert(!_it.valueless_by_exception());
  auto old_it = *this;
  std::visit([&](auto &it) { std::ignore = --it; }, _it);
  return old_it;
}

auto BinTable::iterator::operator-=(difference_type i) -> iterator & {
  assert(!_it.valueless_by_exception());
  std::visit([&](auto &it) { std::ignore = it -= i; }, _it);
  return *this;
}

auto BinTable::iterator::operator-(difference_type i) const -> iterator {
  assert(!_it.valueless_by_exception());
  auto it = *this;
  it -= i;
  return it;
}

auto BinTable::iterator::operator-(const iterator &other) const -> difference_type {
  assert(!_it.valueless_by_exception());
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 - it2;
      },
      _it);
}
}  // namespace hictk
