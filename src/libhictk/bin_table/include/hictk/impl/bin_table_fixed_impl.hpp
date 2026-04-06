// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "hictk/reference.hpp"

namespace hictk {  // NOLINT

template <typename ChromIt>
BinTableFixed::BinTableFixed(ChromIt first_chrom, ChromIt last_chrom, std::uint32_t bin_size,
                             std::size_t bin_offset)
    : BinTableFixed(Reference(first_chrom, last_chrom), bin_size, bin_offset) {}

template <typename ChromNameIt, typename ChromSizeIt>
BinTableFixed::BinTableFixed(ChromNameIt first_chrom_name, ChromNameIt last_chrom_name,
                             ChromSizeIt first_chrom_size, std::uint32_t bin_size,
                             std::size_t bin_offset)
    : BinTableFixed(Reference(first_chrom_name, last_chrom_name, first_chrom_size), bin_size,
                    bin_offset) {}

constexpr std::uint32_t BinTableFixed::resolution() const noexcept { return _bin_size; }

constexpr const Reference &BinTableFixed::chromosomes() const noexcept { return _chroms; }

constexpr const std::vector<std::uint64_t> &BinTableFixed::num_bin_prefix_sum() const noexcept {
  return _num_bins_prefix_sum;
}

constexpr bool BinTableFixed::iterator::operator==(const iterator &other) const noexcept {
  // clang-format off
  return _bin_table == other._bin_table &&
         _chrom_id == other._chrom_id &&
         _rel_bin_id == other._rel_bin_id;
  // clang-format on
}

constexpr bool BinTableFixed::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

constexpr bool BinTableFixed::iterator::operator<(const iterator &other) const noexcept {
  return bin_id() < other.bin_id();
}

constexpr bool BinTableFixed::iterator::operator<=(const iterator &other) const noexcept {
  return bin_id() <= other.bin_id();
}

constexpr bool BinTableFixed::iterator::operator>(const iterator &other) const noexcept {
  return bin_id() > other.bin_id();
}

constexpr bool BinTableFixed::iterator::operator>=(const iterator &other) const noexcept {
  return bin_id() >= other.bin_id();
}

constexpr std::uint32_t BinTableFixed::iterator::resolution() const noexcept {
  assert(_bin_table);
  return _bin_table->resolution();
}

constexpr std::size_t BinTableFixed::iterator::bin_id() const noexcept {
  return _chrom_bin_id + _rel_bin_id;
}

}  // namespace hictk
