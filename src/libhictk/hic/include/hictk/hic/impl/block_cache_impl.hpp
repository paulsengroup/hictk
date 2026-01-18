// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>

#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

constexpr bool BlockID::operator==(const BlockID &other) const noexcept {
  return chrom1_id == other.chrom1_id && chrom2_id == other.chrom2_id && id == other.id;
}

constexpr std::size_t BlockCache::capacity() const noexcept { return _capacity; }
constexpr std::size_t BlockCache::size() const noexcept { return _size; }
constexpr std::size_t BlockCache::capacity_bytes() const noexcept {
  return capacity() * sizeof(ThinPixel<float>);
}
constexpr std::size_t BlockCache::size_bytes() const noexcept {
  return size() * sizeof(ThinPixel<float>);
}

constexpr double BlockCache::hit_rate() const noexcept {
  if (_hits + _misses == 0) {
    return 0.0;
  }
  return static_cast<double>(_hits) / static_cast<double>(_hits + _misses);
}

constexpr void BlockCache::reset_stats() noexcept {
  _hits = 0;
  _misses = 0;
}

constexpr std::size_t BlockCache::hits() const noexcept { return _hits; }
constexpr std::size_t BlockCache::misses() const noexcept { return _misses; }

}  // namespace hictk::hic::internal
