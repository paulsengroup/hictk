// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>

#include "hictk/bin_table_fixed.hpp"
#include "hictk/bin_table_variable.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/reference.hpp"
namespace hictk {

class BinTable {
  std::variant<BinTableFixed, BinTableVariable<>> _table{BinTableFixed{}};

 public:
  enum class Type : std::uint_fast8_t { fixed, variable };
  using BinTableVar = decltype(_table);
  class iterator;
  friend iterator;

  BinTable() = default;
  template <typename BinTableT>
  explicit BinTable(BinTableT table);
  BinTable(Reference chroms, std::uint32_t bin_size, std::size_t bin_offset = 0);
  template <typename ChromIt>
  BinTable(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size,
           std::size_t bin_offset = 0);
  template <typename ChromNameIt, typename ChromSizeIt>
  BinTable(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name, ChromSizeIt first_chrom_size,
           std::uint32_t bin_size, std::size_t bin_offset = 0);
  template <typename I>
  BinTable(Reference chroms, const std::vector<I> &start_pos, const std::vector<I> &end_pos,
           I bin_offset = 0);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t num_chromosomes() const;
  [[nodiscard]] constexpr std::uint32_t resolution() const noexcept;
  [[nodiscard]] constexpr const Reference &chromosomes() const noexcept;
  // clang-format off
  [[deprecated("superseded by BinTable::type() == BinTable::Type::fixed")]]
  [[nodiscard]] constexpr bool has_fixed_resolution() const noexcept;
  // clang-format on
  [[nodiscard]] constexpr auto type() const noexcept -> Type;

  [[nodiscard]] constexpr const std::vector<std::uint64_t> &num_bin_prefix_sum() const noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;
  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  [[nodiscard]] BinTable subset(const Chromosome &chrom) const;
  [[nodiscard]] BinTable subset(std::string_view chrom_name) const;
  [[nodiscard]] BinTable subset(std::uint32_t chrom_id) const;

  [[nodiscard]] auto find_overlap(const GenomicInterval &query) const
      -> std::pair<BinTable::iterator, BinTable::iterator>;
  [[nodiscard]] auto find_overlap(const Chromosome &chrom, std::uint32_t start, std::uint32_t end)
      const -> std::pair<BinTable::iterator, BinTable::iterator>;
  [[nodiscard]] auto find_overlap(std::string_view chrom_name, std::uint32_t start,
                                  std::uint32_t end) const
      -> std::pair<BinTable::iterator, BinTable::iterator>;
  [[nodiscard]] auto find_overlap(std::uint32_t chrom_id, std::uint32_t start, std::uint32_t end)
      const -> std::pair<BinTable::iterator, BinTable::iterator>;

  // Map bin_id to Bin
  [[nodiscard]] Bin at(std::uint64_t bin_id) const;
  [[nodiscard]] std::pair<Bin, Bin> at(const GenomicInterval &gi) const;
  [[nodiscard]] Bin at(const Chromosome &chrom, std::uint32_t pos = 0) const;
  [[nodiscard]] Bin at(std::string_view chrom_name, std::uint32_t pos = 0) const;
  [[nodiscard]] Bin at(std::uint32_t chrom_id, std::uint32_t pos) const;
  [[nodiscard]] Bin at_hint(std::uint64_t bin_id, const Chromosome &chrom) const;

  // Map genomic coords to bin_id
  [[nodiscard]] std::pair<std::uint64_t, std::uint64_t> map_to_bin_ids(
      const GenomicInterval &gi) const;
  [[nodiscard]] std::uint64_t map_to_bin_id(const Chromosome &chrom, std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t map_to_bin_id(std::string_view chrom_name, std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t map_to_bin_id(std::uint32_t chrom_id, std::uint32_t pos) const;

  [[nodiscard]] bool operator==(const BinTable &other) const;
  [[nodiscard]] bool operator!=(const BinTable &other) const;

  template <typename BinTableT>
  [[nodiscard]] constexpr const BinTableT &get() const;
  template <typename BinTableT>
  [[nodiscard]] constexpr BinTableT &get();
  [[nodiscard]] constexpr auto get() const noexcept -> const BinTableVar &;
  [[nodiscard]] constexpr auto get() noexcept -> BinTableVar &;

  class iterator {
    friend BinTable;
    std::variant<BinTableFixed::iterator, BinTableVariable<>::iterator> _it{
        BinTableFixed::iterator{}};

    template <typename It>
    explicit iterator(It it) noexcept;

   public:
    using IteratorVar = decltype(_it);
    using difference_type = std::ptrdiff_t;
    using value_type = Bin;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::random_access_iterator_tag;

    constexpr iterator() noexcept = default;

    constexpr bool operator==(const iterator &other) const;
    constexpr bool operator!=(const iterator &other) const;
    constexpr bool operator<(const iterator &other) const;
    constexpr bool operator<=(const iterator &other) const;
    constexpr bool operator>(const iterator &other) const;
    constexpr bool operator>=(const iterator &other) const;

    auto operator*() const -> value_type;
    auto operator[](difference_type i) const -> value_type;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
    auto operator+=(std::size_t i) -> iterator &;
    auto operator+(std::size_t i) const -> iterator;

    auto operator--() -> iterator &;
    auto operator--(int) -> iterator;
    auto operator-=(std::size_t i) -> iterator &;
    auto operator-(std::size_t i) const -> iterator;
    auto operator-(const iterator &other) const -> difference_type;

    template <typename IteratorT>
    [[nodiscard]] constexpr const IteratorT &get() const;
    template <typename IteratorT>
    [[nodiscard]] constexpr IteratorT &get();
    [[nodiscard]] constexpr auto get() const noexcept -> const IteratorVar &;
    [[nodiscard]] constexpr auto get() noexcept -> IteratorVar &;
  };
};

}  // namespace hictk

#include "./impl/bin_table_impl.hpp"  // NOLINT
