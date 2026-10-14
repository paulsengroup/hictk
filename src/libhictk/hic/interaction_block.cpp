// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/interaction_block.hpp"

#include <cstddef>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

InteractionBlock::InteractionBlock(std::size_t id_, [[maybe_unused]] std::size_t block_bin_count,
                                   std::vector<ThinPixel<float>> pixels)
    : _id(id_), _interactions(std::move(pixels)) {}

auto InteractionBlock::operator()() const noexcept -> const BuffT& { return _interactions; }

auto InteractionBlock::begin() const noexcept -> const_iterator { return _interactions.begin(); }
auto InteractionBlock::end() const noexcept -> const_iterator { return _interactions.end(); }
auto InteractionBlock::cbegin() const noexcept -> const_iterator { return begin(); }
auto InteractionBlock::cend() const noexcept -> const_iterator { return end(); }

std::size_t InteractionBlock::id() const noexcept { return _id; }

std::size_t InteractionBlock::size() const noexcept { return _interactions.size(); }

}  // namespace hictk::hic::internal
