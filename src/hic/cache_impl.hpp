// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <numeric>
#include <vector>

namespace hictk::internal {

inline InteractionBlock::InteractionBlock(std::vector<SerializedPixel> interactions) noexcept(
    ndebug_defined()) {
  if (interactions.empty()) {
    return;
  }

  for (auto &&interaction : interactions) {
    _first_col = (std::min)(_first_col, interaction.bin2_id);
    _last_col = (std::max)(_last_col, interaction.bin2_id);
    auto [node, inserted] = this->_interactions.try_emplace(
        interaction.bin1_id, std::vector<SerializedPixel>{std::move(interaction)});
    if (!inserted) {
      node->second.emplace_back(std::move(interaction));
    }
  }
  if constexpr (ndebug_not_defined()) {
    for (auto &[_, buff] : this->_interactions) {
      if (!std::is_sorted(buff.begin(), buff.end())) {
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

inline auto InteractionBlock::find_overlap(std::int64_t first_row,
                                           std::int64_t last_row) const noexcept
    -> std::pair<const_iterator, const_iterator> {
  assert(first_row <= last_row);
  return std::make_pair(_interactions.lower_bound(first_row), _interactions.upper_bound(last_row));
}

inline std::size_t InteractionBlock::size() const noexcept { return _interactions.size(); }

inline std::size_t InteractionBlock::size_in_bytes() const noexcept {
  return sizeof(Pixel<float>) * size();
}

inline std::int64_t InteractionBlock::first_row() const noexcept {
  if (_interactions.empty()) {
    return 0;
  }
  return _interactions.begin()->second.front().bin1_id;
}
inline std::int64_t InteractionBlock::last_row() const noexcept {
  if (_interactions.empty()) {
    return 0;
  }
  return (--_interactions.end())->second.back().bin1_id;
}
inline std::int64_t InteractionBlock::first_col() const noexcept { return _first_col; }
inline std::int64_t InteractionBlock::last_col() const noexcept { return _last_col; }

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

}  // namespace hictk::internal
