// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
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
}  // namespace hictk::hic::internal
