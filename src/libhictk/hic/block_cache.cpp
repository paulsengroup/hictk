// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "hictk/hic/cache.hpp"
// clang-format on

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "hictk/hash.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

BlockCache::BlockCache(std::size_t capacity_bytes)
    : _capacity(capacity_bytes / sizeof(ThinPixel<float>)) {}

auto BlockCache::find(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id) -> Value {
  auto match = _map.find({chrom1_id, chrom2_id, block_id});
  if (match != _map.end()) {
    ++_hits;
    return match->second;
  }

  ++_misses;
  return nullptr;
}

auto BlockCache::emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
                         Value block) -> Value {
  while (size() + block->size() > capacity() && !_map.empty()) {
    pop_oldest();
  }

  const BlockID key{chrom1_id, chrom2_id, block_id};
  _queue.push(key);
  _map.emplace(key, block);
  _size += block->size();
  return block;
}

auto BlockCache::emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
                         InteractionBlock &&block) -> Value {
  return emplace(chrom1_id, chrom2_id, block_id,
                 std::make_shared<InteractionBlock>(std::move(block)));
}

bool BlockCache::try_erase(const BlockID &key) {
  auto it = _map.find(key);
  if (it != _map.end()) {
    _size -= it->second->size();
    _map.erase(it);
    return true;
  }
  return false;
}

bool BlockCache::try_erase(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id) {
  return try_erase({chrom1_id, chrom2_id, block_id});
}

void BlockCache::clear() noexcept {
  reset_stats();
  _map.clear();
  while (!_queue.empty()) {
    _queue.pop();
  }
}

std::size_t BlockCache::num_blocks() const noexcept { return _map.size(); }

void BlockCache::set_capacity(std::size_t new_capacity, bool shrink_to_fit) {
  if (shrink_to_fit) {
    while (new_capacity < size_bytes() && size() != 0) {
      pop_oldest();
    }
  }
  _capacity = new_capacity / sizeof(ThinPixel<float>);
}

void BlockCache::pop_oldest() {
  while (!_map.empty()) {
    const auto erased = try_erase(_queue.front());
    _queue.pop();
    if (erased) {
      break;
    }
  }
}

}  // namespace hictk::hic::internal

std::size_t std::hash<hictk::hic::internal::BlockID>::operator()(
    const hictk::hic::internal::BlockID &bid) const noexcept {
  return hictk::internal::hash_combine(0, bid.chrom1_id, bid.chrom2_id, bid.id);
}
