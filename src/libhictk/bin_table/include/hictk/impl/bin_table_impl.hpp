// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/bin_table_fixed.hpp"
#include "hictk/bin_table_variable.hpp"
#include "hictk/type_traits.hpp"

namespace hictk {
template <typename BinTableT>  // NOLINTNEXTLINE(*-avoid-non-const-global-variables)
BinTable::BinTable(BinTableT table) : _table(std::move(table)) {}

template <typename ChromIt>
BinTable::BinTable(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size,
                   std::size_t bin_offset)
    : BinTable(BinTableFixed(first_chrom, last_chrom, bin_size, bin_offset)) {}

template <typename ChromNameIt, typename ChromSizeIt>
BinTable::BinTable(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name,
                   ChromSizeIt first_chrom_size, std::uint32_t bin_size, std::size_t bin_offset)
    : BinTable(BinTableFixed(first_chrom_name, last_chrom_name, first_chrom_size, bin_size,
                             bin_offset)) {}

template <typename I>
BinTable::BinTable(Reference chroms, const std::vector<I> &start_pos, const std::vector<I> &end_pos,
                   I bin_offset)
    : BinTable(BinTableVariable(std::move(chroms), start_pos, end_pos, bin_offset)) {}

// NOLINTNEXTLINE(bugprone-exception-escape)
constexpr std::uint32_t BinTable::resolution() const noexcept {
  assert(!_table.valueless_by_exception());
  if (type() == Type::fixed) {
    return std::get<BinTableFixed>(_table).resolution();
  }
  return 0;
}

// NOLINTNEXTLINE(bugprone-exception-escape)
constexpr const Reference &BinTable::chromosomes() const noexcept {
  assert(!_table.valueless_by_exception());
  return std::visit([&](const auto &t) -> const Reference & { return t.chromosomes(); }, _table);
}

constexpr bool BinTable::has_fixed_resolution() const noexcept { return type() == Type::fixed; }

constexpr auto BinTable::type() const noexcept -> Type {
  assert(!_table.valueless_by_exception());
  return std::holds_alternative<BinTableFixed>(_table) ? Type::fixed : Type::variable;
}

// NOLINTNEXTLINE(bugprone-exception-escape)
constexpr const std::vector<std::uint64_t> &BinTable::num_bin_prefix_sum() const noexcept {
  assert(!_table.valueless_by_exception());
  return std::visit(
      [&](const auto &t) -> const std::vector<std::uint64_t> & { return t.num_bin_prefix_sum(); },
      _table);
}

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
BinTable::iterator::iterator(It it) noexcept : _it(it) {}

constexpr bool BinTable::iterator::operator==(const iterator &other) const {
  assert(!_it.valueless_by_exception());
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
  assert(!_it.valueless_by_exception());
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 < it2;
      },
      _it);
}

constexpr bool BinTable::iterator::operator<=(const iterator &other) const {
  assert(!_it.valueless_by_exception());
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 <= it2;
      },
      _it);
}

constexpr bool BinTable::iterator::operator>(const iterator &other) const {
  assert(!_it.valueless_by_exception());
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 > it2;
      },
      _it);
}

constexpr bool BinTable::iterator::operator>=(const iterator &other) const {
  assert(!_it.valueless_by_exception());
  return std::visit(
      [&](const auto &it1) {
        const auto &it2 = std::get<remove_cvref_t<decltype(it1)>>(other._it);
        return it1 >= it2;
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
