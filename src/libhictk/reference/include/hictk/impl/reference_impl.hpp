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
    : _buff(construct_chrom_buffer(first_chrom, last_chrom)),
      _map(construct_chrom_map(_buff)),
      _size_prefix_sum(compute_size_prefix_sum(_buff)),
      _longest_chrom(find_longest_chromosome(_buff)),
      _chrom_with_longest_name(find_chromosome_with_longest_name(_buff)) {
  validate();
}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
Reference::Reference(ChromosomeNameIt first_chrom_name, ChromosomeNameIt last_chrom_name,
                     ChromosomeSizeIt first_chrom_size)
    : _buff(construct_chrom_buffer(first_chrom_name, last_chrom_name, first_chrom_size)),
      _map(construct_chrom_map(_buff)),
      _size_prefix_sum(compute_size_prefix_sum(_buff)),
      _longest_chrom(find_longest_chromosome(_buff)),
      _chrom_with_longest_name(find_chromosome_with_longest_name(_buff)) {
  validate();
}

constexpr const std::vector<std::uint64_t>& Reference::chrom_size_prefix_sum() const noexcept {
  return _size_prefix_sum;
}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
auto Reference::construct_chrom_buffer(ChromosomeNameIt first_chrom_name,
                                       ChromosomeNameIt last_chrom_name,
                                       ChromosomeSizeIt first_chrom_size) -> ChromBuff {
  ChromBuff buff{};
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
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  std::for_each(first_chrom, last_chrom, [&](const Chromosome& chrom) {
    chrom_names.emplace_back(chrom.name());
    chrom_sizes.emplace_back(chrom.size());
  });

  return construct_chrom_buffer(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());
}

}  // namespace hictk
