// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <variant>
#include <vector>

#include "hictk/bin.hpp"
#include "hictk/common.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/reference.hpp"
#include "hictk/type_traits.hpp"

namespace hictk {
template <typename BinTableT>  // NOLINTNEXTLINE(*-avoid-non-const-global-variables)
inline BinTable::BinTable(BinTableT table) : _table(std::move(table)) {}

inline BinTable::BinTable(Reference chroms, std::uint32_t bin_size, std::size_t bin_offset)
    : BinTable(BinTableFixed(std::move(chroms), bin_size, bin_offset)) {}

template <typename ChromIt>
inline BinTable::BinTable(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size,
                          std::size_t bin_offset)
    : BinTable(BinTableFixed(first_chrom, last_chrom, bin_size, bin_offset)) {}

template <typename ChromNameIt, typename ChromSizeIt>
inline BinTable::BinTable(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name,
                          ChromSizeIt first_chrom_size, std::uint32_t bin_size,
                          std::size_t bin_offset)
    : BinTable(BinTableFixed(first_chrom_name, last_chrom_name, first_chrom_size, bin_size,
                             bin_offset)) {}

template <typename I>
inline BinTable::BinTable(Reference chroms, const std::vector<I> &start_pos,
                          const std::vector<I> &end_pos, I bin_offset)
    : BinTable(BinTableVariable(std::move(chroms), start_pos, end_pos, bin_offset)) {}

inline std::size_t BinTable::size() const noexcept {
  return std::visit([&](const auto &t) { return t.size(); }, _table);
}

inline bool BinTable::empty() const noexcept { return size() == 0; }

inline std::size_t BinTable::num_chromosomes() const {
  return std::visit([&](const auto &t) { return t.num_chromosomes(); }, _table);
}

constexpr std::uint32_t BinTable::resolution() const noexcept {
  if (type() == BinTable::Type::fixed) {
    return std::get<BinTableFixed>(_table).resolution();
  }
  return 0;
}

constexpr const Reference &BinTable::chromosomes() const noexcept {
  return std::visit([&](const auto &t) -> const Reference & { return t.chromosomes(); }, _table);
}

constexpr bool BinTable::has_fixed_resolution() const noexcept { return type() == Type::fixed; }

constexpr auto BinTable::type() const noexcept -> Type {
  return std::holds_alternative<BinTableFixed>(_table) ? Type::fixed : Type::variable;
}

constexpr const std::vector<std::uint64_t> &BinTable::num_bin_prefix_sum() const noexcept {
  return std::visit(
      [&](const auto &t) -> const std::vector<std::uint64_t> & { return t.num_bin_prefix_sum(); },
      _table);
}

inline auto BinTable::begin() const -> iterator {
  return std::visit([&](const auto &t) { return iterator{t.begin()}; }, _table);
}

inline auto BinTable::end() const -> iterator {
  return std::visit([&](const auto &t) { return iterator{t.end()}; }, _table);
}

inline auto BinTable::cbegin() const -> iterator { return begin(); }

inline auto BinTable::cend() const -> iterator { return end(); }

inline BinTable BinTable::subset(const Chromosome &chrom) const {
  return std::visit([&](const auto &t) { return BinTable{t.subset(chrom)}; }, _table);
}

inline BinTable BinTable::subset(std::string_view chrom_name) const {
  return subset(chromosomes().at(chrom_name));
}

inline BinTable BinTable::subset(std::uint32_t chrom_id) const {
  return subset(chromosomes().at(chrom_id));
}

inline auto BinTable::find_overlap(const GenomicInterval &query) const
    -> std::pair<BinTable::iterator, BinTable::iterator> {
  return std::visit(
      [&](const auto &t) {
        auto its = t.find_overlap(query);
        DISABLE_WARNING_PUSH
        DISABLE_WARNING_MAYBE_UNINITIALIZED
        return std::make_pair(iterator{its.first}, iterator{its.second});
        DISABLE_WARNING_POP
      },
      _table);
}

inline auto BinTable::find_overlap(const Chromosome &chrom, std::uint32_t start, std::uint32_t end)
    const -> std::pair<BinTable::iterator, BinTable::iterator> {
  return find_overlap(GenomicInterval{chrom, start, end});
}

inline auto BinTable::find_overlap(std::string_view chrom_name, std::uint32_t start,
                                   std::uint32_t end) const
    -> std::pair<BinTable::iterator, BinTable::iterator> {
  return find_overlap(chromosomes().at(chrom_name), start, end);
}

inline auto BinTable::find_overlap(std::uint32_t chrom_id, std::uint32_t start, std::uint32_t end)
    const -> std::pair<BinTable::iterator, BinTable::iterator> {
  return find_overlap(chromosomes().at(chrom_id), start, end);
}

// Map bin_id to Bin
inline Bin BinTable::at(std::uint64_t bin_id) const {
  return std::visit([&](const auto &t) { return t.at(bin_id); }, _table);
}

inline std::pair<Bin, Bin> BinTable::at(const GenomicInterval &gi) const {
  return std::visit([&](const auto &t) { return t.at(gi); }, _table);
}

inline Bin BinTable::at(const Chromosome &chrom, std::uint32_t pos) const {
  return std::visit([&](const auto &t) { return t.at(chrom, pos); }, _table);
}

inline Bin BinTable::at(std::string_view chrom_name, std::uint32_t pos) const {
  return std::visit([&](const auto &t) { return t.at(chrom_name, pos); }, _table);
}

inline Bin BinTable::at(std::uint32_t chrom_id, std::uint32_t pos) const {
  return std::visit([&](const auto &t) { return t.at(chrom_id, pos); }, _table);
}

inline Bin BinTable::at_hint(std::uint64_t bin_id, const Chromosome &chrom) const {
  return std::visit([&](const auto &t) { return t.at_hint(bin_id, chrom); }, _table);
}

// Map genomic coords to bin_id
inline std::pair<std::uint64_t, std::uint64_t> BinTable::map_to_bin_ids(
    const GenomicInterval &gi) const {
  return std::visit([&](const auto &t) { return t.map_to_bin_ids(gi); }, _table);
}

inline std::uint64_t BinTable::map_to_bin_id(const Chromosome &chrom, std::uint32_t pos) const {
  return std::visit([&](const auto &t) { return t.map_to_bin_id(chrom, pos); }, _table);
}

inline std::uint64_t BinTable::map_to_bin_id(std::string_view chrom_name, std::uint32_t pos) const {
  return std::visit([&](const auto &t) { return t.map_to_bin_id(chrom_name, pos); }, _table);
}

inline std::uint64_t BinTable::map_to_bin_id(std::uint32_t chrom_id, std::uint32_t pos) const {
  return std::visit([&](const auto &t) { return t.map_to_bin_id(chrom_id, pos); }, _table);
}

inline bool BinTable::operator==(const BinTable &other) const {
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

inline bool BinTable::operator!=(const BinTable &other) const { return !(*this == other); }

template <typename BinTableT>
constexpr const BinTableT &BinTable::get() const {
  return std::get<BinTableT>(get());
}

template <typename BinTableT>
constexpr BinTableT &BinTable::get() {
  return std::get<BinTableT>(get());
}

constexpr auto BinTable::get() const noexcept -> const BinTableVar & { return _table; }

constexpr auto BinTable::get() noexcept -> BinTableVar & { return _table; }

template <typename It>
inline BinTable::iterator::iterator(It it) noexcept : _it(it) {}

constexpr bool BinTable::iterator::operator==(const iterator &other) const {
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 == it2;
      },
      _it);
}

constexpr bool BinTable::iterator::operator!=(const iterator &other) const {
  return !(*this == other);
}

constexpr bool BinTable::iterator::operator<(const iterator &other) const {
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 < it2;
      },
      _it);
}

constexpr bool BinTable::iterator::operator<=(const iterator &other) const {
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 <= it2;
      },
      _it);
}

constexpr bool BinTable::iterator::operator>(const iterator &other) const {
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 > it2;
      },
      _it);
}

constexpr bool BinTable::iterator::operator>=(const iterator &other) const {
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 >= it2;
      },
      _it);
}

inline auto BinTable::iterator::operator*() const -> value_type {
  return std::visit([&](const auto &it) { return *it; }, _it);
}
inline auto BinTable::iterator::operator[](std::size_t i) const -> iterator {
  std::visit([&](const auto &it) { std::ignore = it + i; }, _it);
  return *this;
}

inline auto BinTable::iterator::operator++() -> iterator & {
  std::visit([&](auto &it) { std::ignore = ++it; }, _it);
  return *this;
}

inline auto BinTable::iterator::operator++(int) -> iterator {
  auto old_it = *this;
  std::visit([&](auto &it) { std::ignore = ++it; }, _it);
  return old_it;
}

inline auto BinTable::iterator::operator+=(std::size_t i) -> iterator & {
  std::visit([&](auto &it) { std::ignore = it += i; }, _it);
  return *this;
}

inline auto BinTable::iterator::operator+(std::size_t i) const -> iterator {
  auto it = *this;
  it += i;
  return it;
}

inline auto BinTable::iterator::operator--() -> iterator & {
  std::visit([&](auto &it) { std::ignore = --it; }, _it);
  return *this;
}

inline auto BinTable::iterator::operator--(int) -> iterator {
  auto old_it = *this;
  std::visit([&](auto &it) { std::ignore = --it; }, _it);
  return old_it;
}

inline auto BinTable::iterator::operator-=(std::size_t i) -> iterator & {
  std::visit([&](auto &it) { std::ignore = it -= i; }, _it);
  return *this;
}

inline auto BinTable::iterator::operator-(std::size_t i) const -> iterator {
  auto it = *this;
  it -= i;
  return it;
}

inline auto BinTable::iterator::operator-(const iterator &other) const -> difference_type {
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 - it2;
      },
      _it);
}

template <typename IteratorT>
constexpr const IteratorT &BinTable::iterator::get() const {
  return std::get<IteratorT>(get());
}

template <typename IteratorT>
constexpr IteratorT &BinTable::iterator::get() {
  return std::get<IteratorT>(get());
}

constexpr auto BinTable::iterator::get() const noexcept -> const IteratorVar & { return _it; }

constexpr auto BinTable::iterator::get() noexcept -> IteratorVar & { return _it; }
}  // namespace hictk
