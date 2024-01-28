// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

class InteractionBlock {
 public:
  using Row = std::vector<ThinPixel<float>>;

 private:
  using BuffT = std::vector<ThinPixel<float>>;
  std::size_t _id{};
  BuffT _interactions{};

 public:
  using iterator = BuffT::iterator;
  using const_iterator = BuffT::const_iterator;

  InteractionBlock() = default;
  InteractionBlock(std::size_t id_, std::size_t block_bin_count,
                   std::vector<ThinPixel<float>> pixels);

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

}  // namespace hictk::hic::internal

#include "./impl/interaction_block_impl.hpp"  // NOLINT
