// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "hictk/hic/cache.hpp"
// clang-format on

#include <cstddef>
#include <cstdint>
#include <memory>
#include <utility>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/chromosome.hpp"

namespace hictk::hic::internal {

auto WeightCache::get_or_init(std::uint32_t chrom_id, balancing::Method norm) -> Value {
  auto key = std::make_pair(chrom_id, std::move(norm));
  auto it = _weights.find(key);
  if (it != _weights.end()) {
    return it->second;
  }

  return _weights.emplace(std::move(key), std::make_shared<balancing::Weights>()).first->second;
}

auto WeightCache::get_or_init(const Chromosome &chrom, balancing::Method norm) -> Value {
  return get_or_init(chrom.id(), std::move(norm));
}

auto WeightCache::at(std::uint32_t chrom_id, balancing::Method norm) const -> Value {
  return _weights.at(std::make_pair(chrom_id, std::move(norm)));
}

auto WeightCache::at(const Chromosome &chrom, balancing::Method norm) const -> Value {
  return at(chrom.id(), std::move(norm));
}

void WeightCache::clear() noexcept { _weights.clear(); }
std::size_t WeightCache::size() const noexcept { return _weights.size(); }

}  // namespace hictk::hic::internal
