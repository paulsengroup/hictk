// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/chromosome.hpp"

namespace hictk::hic::internal {

inline auto WeightCache::get_or_init(std::uint32_t chrom_id, balancing::Method norm) -> Value {
  auto key = std::make_pair(chrom_id, std::move(norm));
  auto it = _weights.find(key);
  if (it != _weights.end()) {
    return it->second;
  }

  return _weights.emplace(std::move(key), std::make_shared<balancing::Weights>()).first->second;
}

inline auto WeightCache::get_or_init(const Chromosome &chrom, balancing::Method norm) -> Value {
  return get_or_init(chrom.id(), std::move(norm));
}

inline auto WeightCache::at(std::uint32_t chrom_id, balancing::Method norm) const -> Value {
  return _weights.at(std::make_pair(chrom_id, std::move(norm)));
}

inline auto WeightCache::at(const Chromosome &chrom, balancing::Method norm) const -> Value {
  return at(chrom.id(), std::move(norm));
}

inline void WeightCache::clear() noexcept { _weights.clear(); }
inline std::size_t WeightCache::size() const noexcept { return _weights.size(); }

}  // namespace hictk::hic::internal
