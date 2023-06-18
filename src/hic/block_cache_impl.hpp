// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hictk/hash.hpp"
#include "hictk/hic/footer.hpp"

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
                                          const std::vector<SerializedPixel> &pixels)
    : _id(id_), _size(pixels.size()) {
  if (pixels.empty()) {
    return;
  }

  for (const SerializedPixel &p : pixels) {
    const auto b1 = static_cast<std::uint64_t>(p.bin1_id);
    const auto b2 = static_cast<std::uint64_t>(p.bin2_id);

    _first_bin1_id = (std::min)(b1, _first_bin1_id);
    _first_bin2_id = (std::min)(b2, _first_bin2_id);
    _last_bin1_id = (std::max)(b1, _last_bin1_id);
    _last_bin2_id = (std::max)(b2, _last_bin2_id);

    auto [node, _] = this->_interactions.try_emplace(b1, Row{});
    node->second.push_back({b2, p.count});
  }
  if constexpr (ndebug_not_defined()) {
    for (auto &[_, buff] : this->_interactions) {
      if (!std::is_sorted(buff.begin(), buff.end(), [](const ThinPixel &p1, const ThinPixel &p2) {
            return p1.bin2_id < p2.bin2_id;
          })) {
        throw std::runtime_error("InteractionBlock is not sorted!");
      }
    }
  }
}

inline auto InteractionBlock::operator()() const noexcept -> const BuffT & { return _interactions; }

inline auto InteractionBlock::begin() noexcept -> iterator { return _interactions.begin(); }

inline auto InteractionBlock::begin() const noexcept -> const_iterator {
  return _interactions.begin();
}

inline auto InteractionBlock::cbegin() const noexcept -> const_iterator { return begin(); }

inline auto InteractionBlock::end() noexcept -> iterator { return _interactions.end(); }

inline auto InteractionBlock::end() const noexcept -> const_iterator { return _interactions.end(); }

inline auto InteractionBlock::cend() const noexcept -> const_iterator { return end(); }

inline std::size_t InteractionBlock::id() const noexcept { return _id; }
inline const Chromosome &InteractionBlock::chrom1() const noexcept {
  assert(_chrom1);
  return *_chrom1;
}
inline const Chromosome &InteractionBlock::chrom2() const noexcept {
  assert(_chrom2);
  return *_chrom2;
}

inline std::size_t InteractionBlock::first_bin1_id() const noexcept { return _first_bin1_id; }
inline std::size_t InteractionBlock::first_bin2_id() const noexcept { return _first_bin2_id; }
inline std::size_t InteractionBlock::last_bin1_id() const noexcept { return _last_bin1_id; }
inline std::size_t InteractionBlock::last_bin2_id() const noexcept { return _last_bin2_id; }

inline auto InteractionBlock::find(std::uint64_t row) const noexcept -> const_iterator {
  return _interactions.find(row);
}

inline std::size_t InteractionBlock::size() const noexcept { return _size; }

inline std::size_t InteractionBlock::size_in_bytes() const noexcept {
  return sizeof(ThinPixel) * size();
}

constexpr bool BlockID::operator==(const BlockID &other) const noexcept {
  return chrom1_id == other.chrom1_id && chrom2_id == other.chrom2_id && id == other.id;
}

inline BlockCache::BlockCache(std::size_t capacity) : _map(capacity), _capacity(capacity) {
  assert(capacity != 0);
}

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
  while (_size + block->size() > capacity() && !_map.empty()) {
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

constexpr std::size_t BlockCache::capacity() const noexcept { return _capacity; }
constexpr std::size_t BlockCache::size() const noexcept { return _size; }
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

constexpr std::size_t BlockCache::hits() const noexcept { return _hits; }
constexpr std::size_t BlockCache::misses() const noexcept { return _misses; }

inline void BlockCache::pop_oldest() {
  auto it = _map.find(_queue.front());
  _queue.pop();
  _size -= it->second->size();
  _map.erase(it);
}

}  // namespace hictk::hic::internal
