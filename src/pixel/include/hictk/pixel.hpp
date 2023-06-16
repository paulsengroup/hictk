// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <limits>
#include <queue>
#include <string_view>
#include <type_traits>

#include "hictk/bin_table.hpp"

namespace hictk {

class Chromosome;

struct PixelCoordinates {
  Bin bin1;
  Bin bin2;

  PixelCoordinates() = default;
  PixelCoordinates(Bin bin1_, Bin bin2_) noexcept;
  explicit PixelCoordinates(std::pair<Bin, Bin> bins) noexcept;
  explicit PixelCoordinates(Bin bin) noexcept;

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool operator==(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator!=(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator<(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator<=(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator>(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator>=(const PixelCoordinates &other) const noexcept;

  [[nodiscard]] bool is_intra() const noexcept;
};

template <typename N>
struct Pixel {
  static_assert(std::is_arithmetic_v<N>);

  PixelCoordinates coords{};
  N count{};

  Pixel() = default;
  explicit Pixel(Bin bin, N count_ = 0) noexcept;
  Pixel(Bin bin1_, Bin bin2_, N count_ = 0) noexcept;
  explicit Pixel(PixelCoordinates coords_, N count_ = 0) noexcept;
  Pixel(const Chromosome &chrom, std::uint32_t start, std::uint32_t end, N count_ = 0) noexcept;
  Pixel(const Chromosome &chrom1, std::uint32_t start1, std::uint32_t end1,
        const Chromosome &chrom2, std::uint32_t start2, std::uint32_t end2, N count_ = 0) noexcept;

  Pixel(const BinTable &bins, std::uint64_t bin1_id, std::uint64_t bin2_id, N count_ = 0);
  Pixel(const BinTable &bins, std::uint64_t bin_id, N count_ = 0);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool operator==(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator!=(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator<(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator<=(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator>(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator>=(const Pixel<N> &other) const noexcept;
};

namespace internal {

/// This class is basically a wrapper around a priority queue of objects of type Node
/// Node consist of a pixel and an index. The index represent from which iterator the
/// pixel was read. This allows us to know from which iterator we should read the next pixel (i.e.
/// the same iterator from which the top pixel originated)
template <typename PixelIt>
class PixelMerger {
  using N = decltype(std::declval<PixelIt>()->count);
  struct Node {
    Pixel<N> pixel{};  // NOLINT
    std::size_t i{};   // NOLINT

    bool operator<(const Node &other) const noexcept;
    bool operator>(const Node &other) const noexcept;
    bool operator==(const Node &other) const noexcept;
    bool operator!=(const Node &other) const noexcept;
  };

  std::vector<Pixel<N>> _buffer{};
  std::priority_queue<Node, std::vector<Node>, std::greater<>> _pqueue{};

  std::vector<PixelIt> _heads{};
  std::vector<PixelIt> _tails{};

 public:
  PixelMerger() = delete;
  PixelMerger(std::vector<PixelIt> head, std::vector<PixelIt> tail);
  [[nodiscard]] auto next() -> Pixel<N>;

 private:
  void replace_top_node(std::size_t i);
};
}  // namespace internal

}  // namespace hictk

#include "../../pixel_impl.hpp"
