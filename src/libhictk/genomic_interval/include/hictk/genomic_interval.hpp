// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT
#pragma once

#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <string_view>
#include <tuple>

#include "hictk/chromosome.hpp"
#include "hictk/reference.hpp"

namespace hictk {
class GenomicInterval {
  static inline const Chromosome null_chrom{};

  Chromosome _chrom{null_chrom};
  std::uint32_t _start{};
  std::uint32_t _end{};

 public:
  enum class Type : std::uint_fast8_t { BED, UCSC };

  constexpr GenomicInterval() = default;
  explicit GenomicInterval(const Chromosome &chrom_) noexcept;
  GenomicInterval(Chromosome chrom_, std::uint32_t start_, std::uint32_t end) noexcept;
  [[nodiscard]] static std::tuple<std::string, std::uint32_t, std::uint32_t> parse(
      const std::string &query, Type type = Type::UCSC);
  [[nodiscard]] static std::tuple<std::string, std::uint32_t, std::uint32_t> parse_ucsc(
      std::string buffer);
  [[nodiscard]] static std::tuple<std::string, std::uint32_t, std::uint32_t> parse_bed(
      std::string_view buffer, char sep = '\t');
  [[nodiscard]] static GenomicInterval parse(const Reference &chroms, const std::string &query,
                                             Type type = Type::UCSC);
  [[nodiscard]] static GenomicInterval parse_ucsc(const Reference &chroms, std::string buffer);
  [[nodiscard]] static GenomicInterval parse_bed(const Reference &chroms, std::string_view buffer,
                                                 char sep = '\t');

  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] bool operator==(const GenomicInterval &other) const noexcept;
  [[nodiscard]] bool operator!=(const GenomicInterval &other) const noexcept;

  [[nodiscard]] bool operator<(const GenomicInterval &other) const noexcept;
  [[nodiscard]] bool operator<=(const GenomicInterval &other) const noexcept;

  [[nodiscard]] bool operator>(const GenomicInterval &other) const noexcept;
  [[nodiscard]] bool operator>=(const GenomicInterval &other) const noexcept;

  [[nodiscard]] const Chromosome &chrom() const noexcept;
  [[nodiscard]] constexpr std::uint32_t start() const noexcept;
  [[nodiscard]] constexpr std::uint32_t end() const noexcept;
  [[nodiscard]] constexpr std::uint32_t size() const noexcept;
  [[nodiscard]] constexpr bool empty() const noexcept;
};

[[nodiscard]] std::uint64_t area(const GenomicInterval &gi, std::uint32_t resolution = 1,
                                 bool upper_triangular = false) noexcept;
[[nodiscard]] std::uint64_t area(const GenomicInterval &gi1, const GenomicInterval &gi2,
                                 std::uint32_t resolution = 1,
                                 bool upper_triangular = false) noexcept;
[[nodiscard]] std::uint64_t area(std::uint32_t start_pos, std::uint32_t end_pos,
                                 std::uint32_t resolution = 1,
                                 bool upper_triangular = false) noexcept;
[[nodiscard]] std::uint64_t area(std::uint32_t start1_pos, std::uint32_t end1_pos,
                                 std::uint32_t start2_pos, std::uint32_t end2_pos,
                                 std::uint32_t resolution = 1,
                                 bool upper_triangular = false) noexcept;

}  // namespace hictk

namespace std {
template <>
struct hash<hictk::GenomicInterval> {
  std::size_t operator()(const hictk::GenomicInterval &gi) const noexcept;
};
}  // namespace std

#include "./impl/genomic_interval_impl.hpp"  // NOLINT
