// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/string.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/type_traits.hpp"
#include "hictk/weights.hpp"

namespace hictk::cooler {

namespace internal {
[[nodiscard]] inline std::vector<std::uint64_t> import_chrom_offsets(const Dataset &dset,
                                                                     std::size_t expected_size) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto offsets = dset.read_all<std::vector<std::uint64_t>>();
  try {
    if (offsets.size() != expected_size) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("expected {} offsets, found {}"), expected_size, offsets.size()));
    }
    if (offsets.front() != 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("first offset should be 0, found {}"), offsets.front()));
    }
    if (!std::is_sorted(offsets.begin(), offsets.end())) {
      throw std::runtime_error("offsets are not in ascending order");
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to import offsets from {}: {}"), dset.uri(), e.what()));
  }

  return offsets;
}
}  // namespace internal

template <typename ChromIt>
void File::read_index_chunk(ChromIt first_chrom, ChromIt last_chrom) const {
  assert(_index);
  try {
    std::for_each(first_chrom, last_chrom, [&](const Chromosome &chrom) {
      if (_index->size(chrom.id()) != 1) {
        return;
      }

      auto chrom_offset_dset = dataset("indexes/chrom_offset");
      auto bin_offset_dset = dataset("indexes/bin1_offset");
      const auto chrom_offsets =
          internal::import_chrom_offsets(chrom_offset_dset, chromosomes().size() + 1);

      auto offset1 = static_cast<std::ptrdiff_t>(chrom_offsets[chrom.id()]);
      auto offset2 = static_cast<std::ptrdiff_t>(chrom_offsets[chrom.id() + 1]);
      auto first = bin_offset_dset.begin<std::uint64_t>() + offset1;
      auto last = bin_offset_dset.begin<std::uint64_t>() + offset2;
      _index->set(chrom, {first, last});

      try {
        _index->validate(chrom);
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(FMT_STRING("index validation failed: {}"), e.what()));
      }
    });

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to import indexes for cooler at URI: \"{}\": {}"),
                    dataset("indexes/bin1_offset").get_parent().uri(), e.what()));
  }
}

}  // namespace hictk::cooler
