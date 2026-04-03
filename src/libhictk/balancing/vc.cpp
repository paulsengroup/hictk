// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/balancing/vc.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "hictk/balancing/weights.hpp"

namespace hictk::balancing {

balancing::Weights VC::get_weights(bool rescale) const {
  if (!rescale) {
    return {_biases, balancing::Weights::Type::DIVISIVE};
  }

  std::vector<double> biases(_biases.size());
  std::uint64_t chrom_id = 0;
  for (std::size_t i = 0; i < _biases.size(); ++i) {
    if (i >= _chrom_offsets[chrom_id + 1]) {
      chrom_id++;
    }
    biases[i] = _biases[i] * _scale[chrom_id];
  }

  std::transform(biases.begin(), biases.end(), biases.begin(), [](const double n) {
    if (std::isnan(n)) {
      return 1.0;
    }
    return n;
  });

  return {biases, balancing::Weights::Type::DIVISIVE};
}

const std::vector<double>& VC::get_scale() const noexcept { return _scale; }

}  // namespace hictk::balancing
