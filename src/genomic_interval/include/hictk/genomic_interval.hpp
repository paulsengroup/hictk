// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT
#pragma once

#include <cstdint>
#include <functional>

#include "hictk/chromosome.hpp"
#include "hictk/reference.hpp"

namespace hictk {
class GenomicInterval {
  static inline const Chromosome null_chrom{};

  Chromosome _chrom{null_chrom};
  std::uint32_t _start{};
  std::uint32_t _end{};

 public:
  enum class Type { BED, UCSC };

  constexpr GenomicInterval() = default;
  explicit GenomicInterval(const Chromosome &chrom_) noexcept;
  GenomicInterval(const Chromosome &chrom_, std::uint32_t start_, std::uint32_t end) noexcept;
  [[nodiscard]] static GenomicInterval parse(const Reference &chroms, std::string query,
                                             Type type = Type::UCSC);
  [[nodiscard]] static GenomicInterval parse_ucsc(const Reference &chroms, std::string query);
  [[nodiscard]] static GenomicInterval parse_bed(const Reference &chroms, std::string_view query,
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
};
}  // namespace hictk

namespace std {
template <>
struct hash<hictk::GenomicInterval> {
  size_t operator()(const hictk::GenomicInterval &gi) const;
};
}  // namespace std

#include "../../genomic_interval_impl.hpp"
