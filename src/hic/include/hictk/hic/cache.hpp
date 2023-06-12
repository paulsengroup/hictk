// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>
#include <tsl/ordered_map.h>

#include <cstdint>
#include <memory>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"
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

 public:
  using iterator = BuffT::iterator;
  using const_iterator = BuffT::const_iterator;

  struct Overlap {
    const_iterator first{};  // NOLINT
    const_iterator last{};   // NOLINT

    [[nodiscard]] auto begin() const noexcept;
    [[nodiscard]] auto end() const noexcept;
    [[nodiscard]] auto cbegin() const noexcept;
    [[nodiscard]] auto cend() const noexcept;
  };

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

  [[nodiscard]] auto at(std::uint64_t row) const noexcept -> const_iterator;
  [[nodiscard]] auto find_overlap(std::uint64_t first_row, std::uint64_t last_row) const noexcept
      -> Overlap;

  [[nodiscard]] bool has_overlap(std::uint64_t first_row, std::uint64_t last_row) const noexcept;

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::size_t size_in_bytes() const noexcept;
};

struct InteractionBlockCmp {
  using is_transparent = void;

  constexpr bool operator()(const InteractionBlock& a, const InteractionBlock& b) const noexcept;
  constexpr bool operator()(const InteractionBlock& a, std::size_t b_id) const noexcept;
  constexpr bool operator()(std::size_t a_id, const InteractionBlock& b) const noexcept;
};

class BlockLRUCache {
  using MapT = tsl::ordered_map<std::size_t, std::shared_ptr<const InteractionBlock>>;
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

}  // namespace hictk::hic::internal

#include "../../../cache_impl.hpp"
