// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>

#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk {

template <typename N>
inline void ExpectedValuesAggregator::add(const ThinPixel<N> &p) {
  add(Pixel<N>{*_bins, p});
}

template <typename N>
inline void ExpectedValuesAggregator::add(const Pixel<N> &p) {
  const auto count = conditional_static_cast<double>(p.count);
  if (std::isnan(count)) {
    return;
  }

  const auto &chrom1 = p.coords.bin1.chrom();
  const auto &chrom2 = p.coords.bin2.chrom();

  if (p.coords.is_intra()) {
    at(chrom1) += count;
    const auto i = p.coords.bin2.id() - p.coords.bin1.id();
    // skip last bin in chromosome if chromosome size is not a multiple of bin size
    // this is done to mimic HiCTools' behavior
    if (i < _actual_distances.size()) {
      _actual_distances[i] += count;
    }
  } else {
    at(chrom1, chrom2) += count;
  }
}

}  // namespace hictk
