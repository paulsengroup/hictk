// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <utility>

#include "hictk/balancing/weights.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"

namespace hictk::hic::internal {

inline auto WeightCache::find_or_emplace(std::uint32_t chrom_id, balancing::Method norm)
    -> Value {
  auto key = std::make_pair(chrom_id, norm);
  auto it = _weights.find(key);
  if (it != _weights.end()) {
    return it->second;
  }

  return _weights.emplace(key, std::make_shared<balancing::Weights>()).first->second;
}

inline auto WeightCache::find_or_emplace(const Chromosome &chrom, balancing::Method norm)
    -> Value {
  return find_or_emplace(chrom.id(), norm);
}

inline void WeightCache::clear() noexcept { _weights.clear(); }
inline std::size_t WeightCache::size() const noexcept { return _weights.size(); }

}  // namespace hictk::hic::internal
