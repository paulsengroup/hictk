// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

constexpr bool operator<(const InteractionBlock &a, const InteractionBlock &b) noexcept {
  return a < b._id;
}
constexpr bool operator==(const InteractionBlock &a, const InteractionBlock &b) noexcept {
  return a == b._id;
}
constexpr bool operator!=(const InteractionBlock &a, const InteractionBlock &b) noexcept {
  return !(a == b);
}

constexpr bool operator<(const InteractionBlock &a, std::size_t b_id) noexcept {
  return a._id < b_id;
}
constexpr bool operator==(const InteractionBlock &a, std::size_t b_id) noexcept {
  return a._id == b_id;
}
constexpr bool operator!=(const InteractionBlock &a, std::size_t b_id) noexcept {
  return !(a == b_id);
}

constexpr bool operator<(std::size_t a_id, const InteractionBlock &b) noexcept {
  return a_id < b._id;
}
constexpr bool operator==(std::size_t a_id, const InteractionBlock &b) noexcept {
  return a_id == b._id;
}
constexpr bool operator!=(std::size_t a_id, const InteractionBlock &b) noexcept {
  return !(a_id == b);
}

inline InteractionBlock::InteractionBlock(std::size_t id_,
                                          [[maybe_unused]] std::size_t block_bin_count,
                                          std::vector<ThinPixel<float>> pixels)
    : _id(id_), _interactions(std::move(pixels)) {}

inline auto InteractionBlock::operator()() const noexcept -> const BuffT & { return _interactions; }

inline auto InteractionBlock::begin() const noexcept -> const_iterator {
  return _interactions.begin();
}
inline auto InteractionBlock::end() const noexcept -> const_iterator { return _interactions.end(); }
inline auto InteractionBlock::cbegin() const noexcept -> const_iterator { return begin(); }
inline auto InteractionBlock::cend() const noexcept -> const_iterator { return end(); }

inline std::size_t InteractionBlock::id() const noexcept { return _id; }

inline std::size_t InteractionBlock::size() const noexcept { return _interactions.size(); }

constexpr bool BlockID::operator==(const BlockID &other) const noexcept {
  return chrom1_id == other.chrom1_id && chrom2_id == other.chrom2_id && id == other.id;
}

inline BlockCache::BlockCache(std::size_t capacity_bytes)
    : _capacity(capacity_bytes / sizeof(ThinPixel<float>)) {}

inline auto BlockCache::find(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id)
    -> Value {
  auto match = _map.find({chrom1_id, chrom2_id, block_id});
  if (match != _map.end()) {
    ++_hits;
    return match->second;
  }

  ++_misses;
  return nullptr;
}

inline auto BlockCache::emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
                                Value block) -> Value {
  while (size() + block->size() > capacity() && !_map.empty()) {
    pop_oldest();
  }

  BlockID key{chrom1_id, chrom2_id, block_id};
  _queue.push(key);
  _map.emplace(std::move(key), block);
  _size += block->size();
  return block;
}

inline auto BlockCache::emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
                                InteractionBlock &&block) -> Value {
  return emplace(chrom1_id, chrom2_id, block_id,
                 std::make_shared<InteractionBlock>(std::move(block)));
}

inline bool BlockCache::try_erase(const BlockID &key) {
  auto it = _map.find(key);
  if (it != _map.end()) {
    _size -= it->second->size();
    _map.erase(it);
    return true;
  }
  return false;
}

inline bool BlockCache::try_erase(std::size_t chrom1_id, std::size_t chrom2_id,
                                  std::size_t block_id) {
  return try_erase({chrom1_id, chrom2_id, block_id});
}

inline void BlockCache::clear() noexcept {
  reset_stats();
  _map.clear();
  while (!_queue.empty()) {
    _queue.pop();
  }
}

constexpr std::size_t BlockCache::capacity() const noexcept { return _capacity; }
constexpr std::size_t BlockCache::size() const noexcept { return _size; }
constexpr std::size_t BlockCache::capacity_bytes() const noexcept {
  return capacity() * sizeof(ThinPixel<float>);
}
constexpr std::size_t BlockCache::size_bytes() const noexcept {
  return size() * sizeof(ThinPixel<float>);
}
inline std::size_t BlockCache::num_blocks() const noexcept { return _map.size(); }

constexpr double BlockCache::hit_rate() const noexcept {
  if (_hits + _misses == 0) {
    return 0.0;
  }
  return double(_hits) / double(_hits + _misses);
}

constexpr void BlockCache::reset_stats() noexcept {
  _hits = 0;
  _misses = 0;
}

inline void BlockCache::set_capacity(std::size_t new_capacity, bool shrink_to_fit) {
  if (shrink_to_fit) {
    while (new_capacity < size_bytes() && size() != 0) {
      pop_oldest();
    }
  }
  _capacity = new_capacity / sizeof(ThinPixel<float>);
}

constexpr std::size_t BlockCache::hits() const noexcept { return _hits; }
constexpr std::size_t BlockCache::misses() const noexcept { return _misses; }

inline void BlockCache::pop_oldest() {
  while (!_map.empty()) {
    const auto erased = try_erase(_queue.front());
    _queue.pop();
    if (erased) {
      break;
    }
  }
}

}  // namespace hictk::hic::internal
