// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <parallel_hashmap/phmap.h>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <queue>
#include <utility>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/hash.hpp"
#include "hictk/hic/interaction_block.hpp"

namespace hictk::hic::internal {
struct BlockID {
  std::size_t chrom1_id;  // NOLINT
  std::size_t chrom2_id;  // NOLINT
  std::size_t id;         // NOLINT
  [[nodiscard]] constexpr bool operator==(const BlockID& other) const noexcept;
};
}  // namespace hictk::hic::internal

template <>
struct std::hash<hictk::hic::internal::BlockID> {
  inline std::size_t operator()(hictk::hic::internal::BlockID const& bid) const noexcept {
    return hictk::internal::hash_combine(0, bid.chrom1_id, bid.chrom2_id, bid.id);
  }
};

namespace hictk::hic::internal {

class BlockCache {
  using Value = std::shared_ptr<const InteractionBlock>;
  std::queue<BlockID> _queue{};
  phmap::flat_hash_map<BlockID, Value> _map{};

  std::size_t _hits{};
  std::size_t _misses{};

  std::size_t _capacity{};
  std::size_t _size{};

 public:
  BlockCache() = delete;
  explicit BlockCache(std::size_t capacity_bytes);

  [[nodiscard]] auto find(std::size_t chrom1_id, std::size_t chrom2_id,
                          std::size_t block_id) -> Value;

  auto emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
               Value block) -> Value;
  auto emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
               InteractionBlock&& block) -> Value;

  bool try_erase(const BlockID& key);
  bool try_erase(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id);
  void clear() noexcept;

  [[nodiscard]] constexpr std::size_t capacity() const noexcept;
  [[nodiscard]] constexpr std::size_t size() const noexcept;
  [[nodiscard]] constexpr std::size_t capacity_bytes() const noexcept;
  [[nodiscard]] constexpr std::size_t size_bytes() const noexcept;
  [[nodiscard]] std::size_t num_blocks() const noexcept;

  [[nodiscard]] constexpr double hit_rate() const noexcept;
  [[nodiscard]] constexpr std::size_t hits() const noexcept;
  [[nodiscard]] constexpr std::size_t misses() const noexcept;
  constexpr void reset_stats() noexcept;
  void set_capacity(std::size_t new_capacity, bool shrink_to_fit = false);

 private:
  void pop_oldest();
};

class WeightCache {
  using Value = std::shared_ptr<balancing::Weights>;
  phmap::flat_hash_map<std::pair<std::uint32_t, balancing::Method>, Value> _weights{};

 public:
  WeightCache() = default;

  [[nodiscard]] auto get_or_init(std::uint32_t chrom_id, balancing::Method norm) -> Value;
  [[nodiscard]] auto get_or_init(const Chromosome& chrom, balancing::Method norm) -> Value;

  [[nodiscard]] auto at(std::uint32_t chrom_id, balancing::Method norm) const -> Value;
  [[nodiscard]] auto at(const Chromosome& chrom, balancing::Method norm) const -> Value;

  void clear() noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
};

}  // namespace hictk::hic::internal

#include "./impl/block_cache_impl.hpp"   // NOLINT
#include "./impl/weight_cache_impl.hpp"  // NOLINT
