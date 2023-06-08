// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "coolerpp/common.hpp"
#include "coolerpp/genomic_interval.hpp"

namespace coolerpp {

class Bin {
 public:
  static constexpr std::uint64_t null_id{(std::numeric_limits<std::uint64_t>::max)()};

 private:
  std::uint64_t _id{null_id};
  GenomicInterval _interval{};

 public:
  constexpr Bin() = default;
  Bin(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end) noexcept;
  Bin(std::uint64_t id_, const Chromosome &chrom_, std::uint32_t start_,
      std::uint32_t end_) noexcept;
  explicit Bin(GenomicInterval interval) noexcept;
  Bin(std::uint64_t id, GenomicInterval interval) noexcept;

  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] bool operator==(const Bin &other) const noexcept;
  [[nodiscard]] bool operator!=(const Bin &other) const noexcept;

  [[nodiscard]] bool operator<(const Bin &other) const noexcept;
  [[nodiscard]] bool operator<=(const Bin &other) const noexcept;

  [[nodiscard]] bool operator>(const Bin &other) const noexcept;
  [[nodiscard]] bool operator>=(const Bin &other) const noexcept;

  [[nodiscard]] constexpr std::uint64_t id() const noexcept;
  [[nodiscard]] const GenomicInterval &interval() const noexcept;
  [[nodiscard]] const Chromosome &chrom() const noexcept;
  [[nodiscard]] constexpr std::uint32_t start() const noexcept;
  [[nodiscard]] constexpr std::uint32_t end() const noexcept;

  [[nodiscard]] constexpr bool has_null_id() const noexcept;
};

struct BinTableConcrete {
  std::vector<const Chromosome *> chroms{};
  std::vector<std::uint32_t> bin_starts{};
  std::vector<std::uint32_t> bin_ends{};
};

class BinTable {
  ChromosomeSet _chroms{};
  std::vector<std::uint64_t> _num_bins_prefix_sum{};
  std::uint32_t _bin_size{std::numeric_limits<std::uint32_t>::max()};

 public:
  class iterator;
  friend iterator;

  BinTable() = default;
  BinTable(ChromosomeSet chroms, std::uint32_t bin_size);
  template <typename ChromIt>
  BinTable(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size);
  template <typename ChromNameIt, typename ChromSizeIt>
  BinTable(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name, ChromSizeIt first_chrom_size,
           std::uint32_t bin_size);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t num_chromosomes() const;
  [[nodiscard]] constexpr std::uint32_t bin_size() const noexcept;
  [[nodiscard]] constexpr const ChromosomeSet &chromosomes() const noexcept;

  [[nodiscard]] constexpr const std::vector<std::uint64_t> &num_bin_prefix_sum() const noexcept;

  [[nodiscard]] constexpr auto begin() const -> iterator;
  [[nodiscard]] constexpr auto end() const -> iterator;
  [[nodiscard]] constexpr auto cbegin() const -> iterator;
  [[nodiscard]] constexpr auto cend() const -> iterator;

  [[nodiscard]] BinTable subset(const Chromosome &chrom) const;
  [[nodiscard]] BinTable subset(std::string_view chrom_name) const;
  [[nodiscard]] BinTable subset(std::uint32_t chrom_id) const;

  [[nodiscard]] auto find_overlap(const GenomicInterval &query) const
      -> std::pair<BinTable::iterator, BinTable::iterator>;
  [[nodiscard]] auto find_overlap(const Chromosome &chrom, std::uint32_t start,
                                  std::uint32_t end) const
      -> std::pair<BinTable::iterator, BinTable::iterator>;
  [[nodiscard]] auto find_overlap(std::string_view chrom_name, std::uint32_t start,
                                  std::uint32_t end) const
      -> std::pair<BinTable::iterator, BinTable::iterator>;
  [[nodiscard]] auto find_overlap(std::uint32_t chrom_id, std::uint32_t start,
                                  std::uint32_t end) const
      -> std::pair<BinTable::iterator, BinTable::iterator>;

  // Map bin_id to Bin
  [[nodiscard]] Bin at(std::uint64_t bin_id) const;
  [[nodiscard]] std::pair<Bin, Bin> at(const GenomicInterval &gi) const;
  [[nodiscard]] Bin at(const Chromosome &chrom, std::uint32_t pos) const;
  [[nodiscard]] Bin at(std::string_view chrom_name, std::uint32_t pos) const;
  [[nodiscard]] Bin at(std::uint32_t chrom_id, std::uint32_t pos) const;
  [[nodiscard]] Bin at_hint(std::uint64_t bin_id, const Chromosome &chrom) const;

  // Map genomic coords to bin_id
  [[nodiscard]] std::pair<std::uint64_t, std::uint64_t> map_to_bin_ids(
      const GenomicInterval &gi) const;
  [[nodiscard]] std::uint64_t map_to_bin_id(const Chromosome &chrom, std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t map_to_bin_id(std::string_view chrom_name, std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t map_to_bin_id(std::uint32_t chrom_id, std::uint32_t pos) const;

  [[nodiscard]] BinTableConcrete concretize() const;

  [[nodiscard]] bool operator==(const BinTable &other) const;
  [[nodiscard]] bool operator!=(const BinTable &other) const;

 private:
  [[nodiscard]] static std::vector<std::uint64_t> compute_num_bins_prefix_sum(
      const ChromosomeSet &chroms, std::uint32_t bin_size);

 public:
  class iterator {
    friend BinTable;
    const BinTable *_bin_table{};
    std::size_t _idx{0};
    std::uint32_t _chrom_id{0};

    static constexpr auto npos = (std::numeric_limits<std::size_t>::max)();
    static constexpr auto nchrom = (std::numeric_limits<std::uint32_t>::max)();

    constexpr explicit iterator(const BinTable &bin_table) noexcept;

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

    auto operator*() const -> value_type;
    auto operator[](std::size_t i) const -> iterator;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
    auto operator+=(std::size_t i) -> iterator &;
    auto operator+(std::size_t i) const -> iterator;

    auto operator--() -> iterator &;
    auto operator--(int) -> iterator;
    auto operator-=(std::size_t i) -> iterator &;
    auto operator-(std::size_t i) const -> iterator;
    auto operator-(const iterator &other) const -> difference_type;

   private:
    [[nodiscard]] static constexpr auto make_end_iterator(const BinTable &table) noexcept
        -> iterator;
    [[nodiscard]] const Chromosome &chromosome(std::uint32_t chrom_id) const;
    [[nodiscard]] const Chromosome &chromosome() const;
    [[nodiscard]] constexpr std::uint32_t bin_size() const noexcept;
    [[nodiscard]] std::uint64_t compute_num_bins() const noexcept;
    [[nodiscard]] std::size_t num_chromosomes() const noexcept;
  };
};

}  // namespace coolerpp

namespace std {
template <>
struct hash<coolerpp::Bin> {
  size_t operator()(const coolerpp::Bin &b) const;
};
}  // namespace std

namespace fmt {
template <>
struct formatter<coolerpp::Bin> {
  enum Presentation { bed, raw, ucsc };
  Presentation presentation{Presentation::raw};

  constexpr format_parse_context::iterator parse(format_parse_context &ctx);
  format_context::iterator format(const coolerpp::Bin &b, format_context &ctx) const;
};
}  // namespace fmt

#include "../../bin_table_impl.hpp"
