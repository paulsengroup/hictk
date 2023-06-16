// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>
#include <tsl/ordered_map.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

class InteractionBlock {
 public:
  struct ThinPixel {
    std::uint64_t bin2_id{};
    float count{};
  };

  using Row = std::vector<ThinPixel>;

 private:
  using BuffT = phmap::btree_map<std::uint64_t, Row>;
  std::size_t _id{};
  BuffT _interactions{};
  const Chromosome* _chrom1{};
  const Chromosome* _chrom2{};
  std::size_t _size{};

 public:
  using iterator = BuffT::iterator;
  using const_iterator = BuffT::const_iterator;

  InteractionBlock() = default;
  InteractionBlock(std::size_t id_, const std::vector<SerializedPixel>& pixels);

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

  [[nodiscard]] auto begin() noexcept -> iterator;
  [[nodiscard]] auto begin() const noexcept -> const_iterator;
  [[nodiscard]] auto cbegin() const noexcept -> const_iterator;

  [[nodiscard]] auto end() noexcept -> iterator;
  [[nodiscard]] auto end() const noexcept -> const_iterator;
  [[nodiscard]] auto cend() const noexcept -> const_iterator;

  [[nodiscard]] std::size_t id() const noexcept;
  [[nodiscard]] const Chromosome& chrom1() const noexcept;
  [[nodiscard]] const Chromosome& chrom2() const noexcept;

  [[nodiscard]] auto find(std::uint64_t row) const noexcept -> const_iterator;

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::size_t size_in_bytes() const noexcept;
};

class BlockLRUCache {
 public:
  struct Key {
    std::size_t chrom1_id;  // NOLINT
    std::size_t chrom2_id;  // NOLINT
    std::size_t id;         // NOLINT

    [[nodiscard]] constexpr bool operator==(const Key& other) const noexcept;
  };

 private:
  struct KeyHasher {
    [[nodiscard]] std::size_t operator()(const Key& k) const noexcept;
  };

  using MapT =
      tsl::ordered_map<Key, std::shared_ptr<const InteractionBlock>, KeyHasher, std::equal_to<>>;
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

  [[nodiscard]] auto find(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id)
      -> iterator;

  auto emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
               mapped_type&& block) -> std::pair<iterator, bool>;
  auto emplace(std::size_t chrom1_id, std::size_t chrom2_id, std::size_t block_id,
               InteractionBlock&& block) -> std::pair<iterator, bool>;

  [[nodiscard]] constexpr double hit_rate() const noexcept;
  [[nodiscard]] constexpr std::size_t hits() const noexcept;
  [[nodiscard]] constexpr std::size_t misses() const noexcept;

 private:
  void erase(key_t key);
  void erase(iterator it);
};

}  // namespace hictk::hic::internal

#include "../../../block_cache_impl.hpp"
