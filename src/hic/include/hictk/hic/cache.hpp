// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <tsl/ordered_map.h>

#include <cstdint>
#include <memory>
#include <vector>

#include "hictk/hic/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::internal {

class InteractionBlock {
  using BuffT = std::map<std::int64_t, std::vector<SerializedPixel>>;
  BuffT _interactions{};
  std::int64_t _first_col{};
  std::int64_t _last_col{};

 public:
  using iterator = BuffT::iterator;
  using const_iterator = BuffT::const_iterator;

  InteractionBlock() = default;
  explicit InteractionBlock(std::vector<SerializedPixel> interactions) noexcept;

  [[nodiscard]] auto operator()() const noexcept -> const BuffT&;

  [[nodiscard]] auto begin() noexcept -> iterator;
  [[nodiscard]] auto begin() const noexcept -> const_iterator;
  [[nodiscard]] auto cbegin() const noexcept -> const_iterator;

  [[nodiscard]] auto end() noexcept -> iterator;
  [[nodiscard]] auto end() const noexcept -> const_iterator;
  [[nodiscard]] auto cend() const noexcept -> const_iterator;

  [[nodiscard]] auto find_overlap(std::int64_t first_row, std::int64_t last_row) const noexcept
      -> std::pair<const_iterator, const_iterator>;

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::size_t size_in_bytes() const noexcept;

  [[nodiscard]] std::int64_t first_row() const noexcept;
  [[nodiscard]] std::int64_t last_row() const noexcept;
  [[nodiscard]] std::int64_t first_col() const noexcept;
  [[nodiscard]] std::int64_t last_col() const noexcept;
};

class BlockLRUCache {
  using MapT = tsl::ordered_map<std::size_t, std::shared_ptr<InteractionBlock>>;
  using key_t = MapT::key_type;
  using mapped_type = MapT::mapped_type;
  using iterator = MapT::iterator;
  using const_iterator = MapT::const_iterator;
  MapT _cache{};

  std::size_t _hits{};
  std::size_t _misses{};

  std::size_t _current_size_bytes{};
  std::size_t _max_size_bytes{500UL * 1024UL * 1024UL};  // 500MB

 public:
  BlockLRUCache() = default;
  explicit BlockLRUCache(std::size_t max_size_in_bytes);

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] constexpr std::size_t size_in_bytes() const noexcept;
  [[nodiscard]] constexpr std::size_t max_size_in_bytes() const noexcept;
  void reset() noexcept;

  [[nodiscard]] auto begin() noexcept -> iterator;
  [[nodiscard]] auto begin() const noexcept -> const_iterator;
  [[nodiscard]] auto cbegin() const noexcept -> const_iterator;

  [[nodiscard]] auto end() noexcept -> iterator;
  [[nodiscard]] auto end() const noexcept -> const_iterator;
  [[nodiscard]] auto cend() const noexcept -> const_iterator;

  [[nodiscard]] auto find(key_t key) -> iterator;

  auto emplace(key_t key, mapped_type&& block) -> std::pair<iterator, bool>;
  auto emplace(key_t key, InteractionBlock&& block) -> std::pair<iterator, bool>;

  [[nodiscard]] constexpr double hit_rate() const noexcept;
  [[nodiscard]] constexpr std::size_t hits() const noexcept;
  [[nodiscard]] constexpr std::size_t misses() const noexcept;

 private:
  void erase(key_t key);
  void erase(iterator it);
};

}  // namespace hictk::internal

#include "../../../cache_impl.hpp"
