// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <nonstd/span.hpp>
#include <queue>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/pixel.hpp"

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

class InteractionBlock {
 public:
  using Row = std::vector<SerializedPixel>;

 private:
  using BuffT = std::vector<SerializedPixel>;
  std::size_t _id{};
  BuffT _interactions{};

 public:
  using iterator = BuffT::iterator;
  using const_iterator = BuffT::const_iterator;

  InteractionBlock() = default;
  InteractionBlock(std::size_t id_, std::size_t block_bin_count,
                   std::vector<SerializedPixel> pixels);

  friend constexpr bool operator<(const InteractionBlock& a, const InteractionBlock& b) noexcept;
  friend constexpr bool operator==(const InteractionBlock& a, const InteractionBlock& b) noexcept;
  friend constexpr bool operator!=(const InteractionBlock& a, const InteractionBlock& b) noexcept;

  friend constexpr bool operator<(const InteractionBlock& a, std::size_t b_id) noexcept;
  friend constexpr bool operator==(const InteractionBlock& a, std::size_t b_id) noexcept;
  friend constexpr bool operator!=(const InteractionBlock& a, std::size_t b_id) noexcept;

  friend constexpr bool operator<(std::size_t a_id, const InteractionBlock& b) noexcept;
  friend constexpr bool operator==(std::size_t a_id, const InteractionBlock& b) noexcept;
  friend constexpr bool operator!=(std::size_t a_id, const InteractionBlock& b) noexcept;

  [[nodiscard]] auto operator()() const noexcept -> const BuffT&;

  [[nodiscard]] auto begin() const noexcept -> const_iterator;
  [[nodiscard]] auto end() const noexcept -> const_iterator;

  [[nodiscard]] auto cbegin() const noexcept -> const_iterator;
  [[nodiscard]] auto cend() const noexcept -> const_iterator;

  [[nodiscard]] std::size_t id() const noexcept;

  [[nodiscard]] std::size_t size() const noexcept;
};

class BlockCache {
  using Value = std::shared_ptr<const InteractionBlock>;
  using MapT = phmap::flat_hash_map<BlockID, Value>;
  std::queue<BlockID> _queue{};
  phmap::flat_hash_map<BlockID, Value> _map{};

  std::size_t _hits{};
  std::size_t _misses{};

  std::size_t _capacity{};
  std::size_t _size{};

 public:
  BlockCache() = delete;
  explicit BlockCache(std::size_t capacity);

  [[nodiscard]] auto find(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id)
      -> Value;

  auto emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id, Value block)
      -> Value;
  auto emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
               InteractionBlock&& block) -> Value;

  bool try_erase(const BlockID& key);
  bool try_erase(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id);
  void clear() noexcept;

  [[nodiscard]] constexpr std::size_t capacity() const noexcept;
  [[nodiscard]] constexpr std::size_t size() const noexcept;
  [[nodiscard]] std::size_t num_blocks() const noexcept;

  [[nodiscard]] constexpr double hit_rate() const noexcept;
  [[nodiscard]] constexpr std::size_t hits() const noexcept;
  [[nodiscard]] constexpr std::size_t misses() const noexcept;
  constexpr void reset_stats() noexcept;
  void set_capacity(std::size_t new_capacity);

 private:
  void pop_oldest();
};

}  // namespace hictk::hic::internal

#include "../../../block_cache_impl.hpp"
