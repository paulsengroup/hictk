// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "hictk/bin.hpp"
#include "hictk/common.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/reference.hpp"

namespace hictk {

template <typename I = std::uint32_t>
class BinTableVariable {
  Reference _chroms{};
  std::vector<I> _bin_end_prefix_sum{0};
  std::vector<std::uint64_t> _num_bins_prefix_sum{0};

 public:
  class iterator;
  friend iterator;

  BinTableVariable() = default;
  BinTableVariable(Reference chroms, const std::vector<I> &start_pos, const std::vector<I> &end_pos,
                   I bin_offset = 0);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t num_chromosomes() const;
  [[nodiscard]] constexpr const Reference &chromosomes() const noexcept;

  [[nodiscard]] constexpr const std::vector<std::uint64_t> &num_bin_prefix_sum() const noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;
  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  [[nodiscard]] BinTableVariable subset(const Chromosome &chrom) const;
  [[nodiscard]] BinTableVariable subset(std::string_view chrom_name) const;
  [[nodiscard]] BinTableVariable subset(std::uint32_t chrom_id) const;

  [[nodiscard]] auto find_overlap(const GenomicInterval &query) const
      -> std::pair<BinTableVariable::iterator, BinTableVariable::iterator>;
  [[nodiscard]] auto find_overlap(const Chromosome &chrom, std::uint32_t start,
                                  std::uint32_t end) const
      -> std::pair<BinTableVariable::iterator, BinTableVariable::iterator>;
  [[nodiscard]] auto find_overlap(std::string_view chrom_name, std::uint32_t start,
                                  std::uint32_t end) const
      -> std::pair<BinTableVariable::iterator, BinTableVariable::iterator>;
  [[nodiscard]] auto find_overlap(std::uint32_t chrom_id, std::uint32_t start,
                                  std::uint32_t end) const
      -> std::pair<BinTableVariable::iterator, BinTableVariable::iterator>;

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

  [[nodiscard]] bool operator==(const BinTableVariable &other) const;
  [[nodiscard]] bool operator!=(const BinTableVariable &other) const;

 private:
  static void validate_bin_coords(const std::vector<I> &start_pos, const std::vector<I> &end_pos);

 public:
  class iterator {
    friend BinTableVariable;
    Bin _value{};
    const BinTableVariable *_bin_table{};
    std::uint32_t _chrom_id{};
    std::uint64_t _bin_id{};

    static constexpr auto null_bin_id = (std::numeric_limits<std::size_t>::max)();
    static constexpr auto nchrom = (std::numeric_limits<std::uint32_t>::max)();

    explicit iterator(const BinTableVariable &bin_table) noexcept;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Bin;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::random_access_iterator_tag;

    constexpr iterator() noexcept = default;

    constexpr bool operator==(const iterator &other) const noexcept;
    constexpr bool operator!=(const iterator &other) const noexcept;
    constexpr bool operator<(const iterator &other) const noexcept;
    constexpr bool operator<=(const iterator &other) const noexcept;
    constexpr bool operator>(const iterator &other) const noexcept;
    constexpr bool operator>=(const iterator &other) const noexcept;

    auto operator*() const noexcept -> value_type;
    auto operator[](difference_type i) const -> value_type;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
    auto operator+=(difference_type i) -> iterator &;
    auto operator+(difference_type i) const -> iterator;

    auto operator--() -> iterator &;
    auto operator--(int) -> iterator;
    auto operator-=(difference_type i) -> iterator &;
    auto operator-(difference_type i) const -> iterator;
    auto operator-(const iterator &other) const -> difference_type;

   private:
    [[nodiscard]] static auto make_end_iterator(const BinTableVariable &table) noexcept -> iterator;
    [[nodiscard]] const Chromosome &chromosome(std::uint32_t chrom_id) const;
    [[nodiscard]] const Chromosome &chromosome() const;
    [[nodiscard]] Bin get_bin() const;
  };
};

}  // namespace hictk

#include "./impl/bin_table_variable_impl.hpp"
