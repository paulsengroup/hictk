// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <numeric>
#include <vector>

#include "hictk/hic/hic_footer.hpp"

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

constexpr bool InteractionBlockCmp::operator()(const InteractionBlock &a,
                                               const InteractionBlock &b) const noexcept {
  return a < b;
}
constexpr bool InteractionBlockCmp::operator()(const InteractionBlock &a,
                                               std::size_t b_id) const noexcept {
  return a < b_id;
}
constexpr bool InteractionBlockCmp::operator()(std::size_t a_id,
                                               const InteractionBlock &b) const noexcept {
  return a_id < b;
}

inline auto InteractionBlock::Overlap::begin() const noexcept { return first; }
inline auto InteractionBlock::Overlap::end() const noexcept { return last; }
inline auto InteractionBlock::Overlap::cbegin() const noexcept { return begin(); }
inline auto InteractionBlock::Overlap::cend() const noexcept { return end(); }

inline std::size_t InteractionBlock::id() const noexcept { return _id; }
inline const Chromosome &InteractionBlock::chrom1() const noexcept {
  assert(_chrom1);
  return *_chrom1;
}
inline const Chromosome &InteractionBlock::chrom2() const noexcept {
  assert(_chrom2);
  return *_chrom2;
}

inline InteractionBlock::InteractionBlock(std::size_t id_,
                                          const std::vector<SerializedPixel> &pixels)
    : _id(id_) {
  if (pixels.empty()) {
    return;
  }

  for (const SerializedPixel &p : pixels) {
    const auto b1 = static_cast<std::uint64_t>(p.bin1_id);
    const auto b2 = static_cast<std::uint64_t>(p.bin2_id);
    auto [node, inserted] = this->_interactions.try_emplace(b1, Row{{b2, p.count}});
    if (!inserted) {
      node->second.emplace_back(ThinPixel{b2, p.count});
    }
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

inline auto InteractionBlock::find(std::uint64_t row) const noexcept -> const_iterator {
  return _interactions.find(row);
}

inline auto InteractionBlock::find_overlap(std::uint64_t first_row,
                                           std::uint64_t last_row) const noexcept -> Overlap {
  assert(first_row <= last_row);
  return {_interactions.lower_bound(first_row), _interactions.upper_bound(last_row)};
}

inline bool InteractionBlock::has_overlap(std::uint64_t first_row,
                                          std::uint64_t last_row) const noexcept {
  auto overlap = find_overlap(first_row, last_row);

  return overlap.begin() != this->_interactions.end();
}

inline std::size_t InteractionBlock::size() const noexcept { return _interactions.size(); }

inline std::size_t InteractionBlock::size_in_bytes() const noexcept {
  return sizeof(Pixel<float>) * size();
}

inline BlockLRUCache::BlockLRUCache(std::size_t max_size_in_bytes)
    : _max_size_bytes(max_size_in_bytes) {
  if (_max_size_bytes == 0) {
    throw std::runtime_error("Invalid block cache capacity: capacity cannot be 0");
  }
}

constexpr std::size_t BlockLRUCache::size_in_bytes() const noexcept { return _current_size_bytes; }
constexpr std::size_t BlockLRUCache::max_size_in_bytes() const noexcept { return _max_size_bytes; }

inline void BlockLRUCache::reset() noexcept {
  _cache.clear();
  _current_size_bytes = 0;
  _hits = 0;
  _misses = 0;
}

inline auto BlockLRUCache::begin() noexcept -> iterator { return _cache.begin(); }

inline auto BlockLRUCache::begin() const noexcept -> const_iterator { return _cache.begin(); }

inline auto BlockLRUCache::cbegin() const noexcept -> const_iterator { return begin(); }

inline auto BlockLRUCache::end() noexcept -> iterator { return _cache.end(); }

inline auto BlockLRUCache::end() const noexcept -> const_iterator { return _cache.end(); }

inline auto BlockLRUCache::cend() const noexcept -> const_iterator { return end(); }

inline auto BlockLRUCache::find(key_t key) -> iterator {
  auto it = _cache.find(key);
  if (it == end()) {
    _misses++;
    return it;
  }
  _hits++;
  auto block = it->second;
  std::ignore = _cache.erase(it);
  auto node = _cache.emplace(key, std::move(block));
  return node.first;
}

inline void BlockLRUCache::erase(key_t key) {
  auto it = _cache.find(key);
  assert(it != _cache.end());
  erase(it);
}

inline void BlockLRUCache::erase(iterator it) {
  _current_size_bytes -= it->second->size_in_bytes();
  std::ignore = _cache.erase(it);
}

inline auto BlockLRUCache::emplace(key_t key, mapped_type &&block) -> std::pair<iterator, bool> {
  assert(block);

  while (size() != 0 && size_in_bytes() + block->size_in_bytes() > max_size_in_bytes()) {
    erase(begin());
  }

  _current_size_bytes += block->size_in_bytes();
  return _cache.emplace(key, std::move(block));
}

inline auto BlockLRUCache::emplace(key_t key, InteractionBlock &&block)
    -> std::pair<iterator, bool> {
  return emplace(key, std::make_shared<InteractionBlock>(std::move(block)));
}

constexpr double BlockLRUCache::hit_rate() const noexcept {
  if (_hits + _misses == 0) {
    return 0.0;
  }
  return double(_hits) / double(_hits + _misses);
}

constexpr std::size_t BlockLRUCache::hits() const noexcept { return _hits; }

constexpr std::size_t BlockLRUCache::misses() const noexcept { return _misses; }

inline std::size_t BlockLRUCache::size() const noexcept { return _cache.size(); }

}  // namespace hictk::hic::internal
