// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"

namespace hictk {

template <typename ChromosomeIt>
Reference::Reference(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
    : Reference(construct_chrom_buffer(first_chrom, last_chrom)) {}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
Reference::Reference(ChromosomeNameIt first_chrom_name, ChromosomeNameIt last_chrom_name,
                     ChromosomeSizeIt first_chrom_size)
    : Reference(construct_chrom_buffer(first_chrom_name, last_chrom_name, first_chrom_size)) {}

constexpr const std::vector<std::uint64_t>& Reference::chrom_size_prefix_sum() const noexcept {
  return _size_prefix_sum;
}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
auto Reference::construct_chrom_buffer(ChromosomeNameIt first_chrom_name,
                                       ChromosomeNameIt last_chrom_name,
                                       ChromosomeSizeIt first_chrom_size) -> ChromBuff {
  const auto num_chroms =
      static_cast<std::size_t>(std::distance(first_chrom_name, last_chrom_name));
  ChromBuff buff{};
  buff.reserve(num_chroms);

  while (first_chrom_name != last_chrom_name) {
    if (std::string_view{*first_chrom_name}.empty()) {
      throw std::runtime_error("found chromosome with empty name");
    }
    buff.emplace_back(static_cast<std::uint32_t>(buff.size()), std::string{*first_chrom_name},
                      conditional_static_cast<std::uint32_t>(*first_chrom_size));

    ++first_chrom_name;  // NOLINT(*-pointer-arithmetic)
    ++first_chrom_size;  // NOLINT(*-pointer-arithmetic)
  }

  return buff;
}

template <typename ChromosomeIt>
auto Reference::construct_chrom_buffer(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
    -> ChromBuff {
  const auto num_chroms = static_cast<std::size_t>(std::distance(first_chrom, last_chrom));
  ChromBuff buff{};
  buff.reserve(num_chroms);
  std::transform(first_chrom, last_chrom, back_inserter(buff),
                 [chrom_id = std::uint32_t{}](const Chromosome& chrom) mutable {
                   if (chrom.name().empty()) {
                     throw std::runtime_error("found chromosome with empty name");
                   }
                   return Chromosome{chrom_id++, chrom.name_ptr(), chrom.size()};
                 });

  return buff;
}

}  // namespace hictk
