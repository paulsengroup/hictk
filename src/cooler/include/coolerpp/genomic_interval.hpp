// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT
#pragma once

#include <cstdint>
#include <functional>

#include "coolerpp/chromosome.hpp"

namespace coolerpp {
class GenomicInterval {
  static inline const Chromosome null_chrom{};

  std::reference_wrapper<const Chromosome> _chrom{null_chrom};
  std::uint32_t _start{};
  std::uint32_t _end{};

 public:
  constexpr GenomicInterval() = default;
  explicit GenomicInterval(const Chromosome &chrom_) noexcept;
  GenomicInterval(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end) noexcept;
  [[nodiscard]] static GenomicInterval parse_ucsc(const ChromosomeSet &chroms, std::string query);
  [[nodiscard]] static GenomicInterval parse_bed(const ChromosomeSet &chroms,
                                                 std::string_view query, char sep = '\t');

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
};
}  // namespace coolerpp

namespace std {
template <>
struct hash<coolerpp::GenomicInterval> {
  size_t operator()(const coolerpp::GenomicInterval &gi) const;
};
}  // namespace std

namespace fmt {
template <>
struct formatter<coolerpp::GenomicInterval> {
  enum Presentation { bed, ucsc };
  Presentation presentation{Presentation::ucsc};

  constexpr format_parse_context::iterator parse(format_parse_context &ctx);
  format_context::iterator format(const coolerpp::GenomicInterval &gi, format_context &ctx) const;
};
}  // namespace fmt

#include "../../genomic_interval_impl.hpp"
