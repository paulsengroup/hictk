// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/balancing/scale.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/weights.hpp"

namespace hictk::balancing {

std::size_t SCALE::size() const noexcept { return _biases.size(); }

void SCALE::reset_iter() noexcept {
  _iter = 0;
  while (!_error_queue_iter.empty()) {
    _error_queue_iter.pop();
  }
}

Weights SCALE::get_weights(bool rescale) const {
  if (!rescale) {
    return {_biases, Weights::Type::DIVISIVE};
  }

  std::vector<double> biases(_biases.size());
  std::uint64_t chrom_id = 0;
  for (std::size_t i = 0; i < _biases.size(); ++i) {
    if (i >= _chrom_offsets[chrom_id + 1]) {
      chrom_id++;
    }
    biases[i] = _biases[i] * _scale[chrom_id];
  }

  return {biases, Weights::Type::DIVISIVE};
}

const std::vector<double>& SCALE::get_scale() const noexcept { return _scale; }

void SCALE::geometric_mean(const std::vector<double>& v1, const std::vector<double>& v2,
                           std::vector<double>& vout) noexcept {
  assert(v1.size() == v2.size());
  assert(vout.size() == v1.size());

  for (std::size_t i = 0; i < v1.size(); ++i) {
    vout[i] = std::sqrt(v1[i] * v2[i]);
  }
}

std::pair<double, std::uint64_t> SCALE::compute_convergence_error(
    const std::vector<double>& biases, const std::vector<double>& current,
    const std::vector<bool>& bad, double tolerance) noexcept {
  assert(biases.size() == current.size());
  assert(biases.size() == bad.size());

  double error = 0;
  std::uint64_t num_fail = 0;
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (bad[i]) {
      continue;
    }
    const auto rel_err = std::abs((biases[i] - current[i]) / (biases[i] + current[i]));
    error = std::max(rel_err, error);
    num_fail += static_cast<std::uint64_t>(rel_err > tolerance);
  }

  return std::make_pair(error, num_fail);
}

double SCALE::compute_final_error(const internal::VectorOfAtomicDecimals& col,
                                  const std::vector<double>& scale,
                                  const std::vector<double>& target,
                                  const std::vector<bool>& bad) noexcept {
  assert(col.size() == scale.size());
  assert(col.size() == target.size());
  assert(col.size() == bad.size());

  double error = 0.0;

  for (std::size_t i = 0; i < col.size(); ++i) {
    if (bad[i]) {
      continue;
    }
    const auto err1 = std::abs((col[i] * scale[i]) - target[i]);
    error = std::max(error, err1);
  }

  return error;
}

void SCALE::multiply(std::vector<double>& v1, const std::vector<double>& v2) noexcept {
  assert(v1.size() == v2.size());
  for (std::size_t i = 0; i < v1.size(); ++i) {
    v1[i] *= v2[i];
  }
}

VC::Type SCALE::map_type_to_vc(Type type) noexcept {
  switch (type) {
    case Type::cis:
      return VC::Type::cis;
    case Type::trans:
      return VC::Type::trans;
    case Type::gw:
      return VC::Type::gw;
  }

  HICTK_UNREACHABLE_CODE;
}

}  // namespace hictk::balancing
