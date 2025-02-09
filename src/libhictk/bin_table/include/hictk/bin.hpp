// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>

#include "hictk/chromosome.hpp"
#include "hictk/genomic_interval.hpp"

namespace hictk {

class Bin {
 public:
  static constexpr std::uint64_t null_id{(std::numeric_limits<std::uint64_t>::max)()};
  static constexpr std::uint32_t rel_null_id{(std::numeric_limits<std::uint32_t>::max)()};

 private:
  std::uint64_t _id{null_id};
  std::uint32_t _rel_id{rel_null_id};
  GenomicInterval _interval{};

 public:
  constexpr Bin() = default;
  Bin(Chromosome chrom_, std::uint32_t start_, std::uint32_t end) noexcept;
  Bin(std::uint64_t id_, std::uint32_t rel_id_, Chromosome chrom_, std::uint32_t start_,
      std::uint32_t end_) noexcept;
  explicit Bin(GenomicInterval interval) noexcept;
  Bin(std::uint64_t id_, std::uint32_t rel_id_, GenomicInterval interval) noexcept;

  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] bool operator==(const Bin &other) const noexcept;
  [[nodiscard]] bool operator!=(const Bin &other) const noexcept;

  [[nodiscard]] bool operator<(const Bin &other) const noexcept;
  [[nodiscard]] bool operator<=(const Bin &other) const noexcept;

  [[nodiscard]] bool operator>(const Bin &other) const noexcept;
  [[nodiscard]] bool operator>=(const Bin &other) const noexcept;

  [[nodiscard]] constexpr std::uint64_t id() const noexcept;
  [[nodiscard]] constexpr std::uint32_t rel_id() const noexcept;
  [[nodiscard]] const GenomicInterval &interval() const noexcept;
  [[nodiscard]] const Chromosome &chrom() const noexcept;
  [[nodiscard]] constexpr std::uint32_t start() const noexcept;
  [[nodiscard]] constexpr std::uint32_t end() const noexcept;

  [[nodiscard]] constexpr bool has_null_id() const noexcept;
};

}  // namespace hictk

namespace std {
template <>
struct hash<hictk::Bin> {
  size_t operator()(const hictk::Bin &b) const;
};
}  // namespace std

#include "./impl/bin_impl.hpp"
