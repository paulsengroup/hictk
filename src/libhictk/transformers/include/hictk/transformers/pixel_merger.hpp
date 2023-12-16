// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <functional>
#include <iterator>
#include <memory>
#include <queue>
#include <type_traits>
#include <vector>

#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::transformers {

/// This class is basically a wrapper around a priority queue of objects of type Node
/// Node consist of a pixel and an index. The index represent from which iterator the
/// pixel was read. This allows us to know from which iterator we should read the next pixel (i.e.
/// the same iterator from which the top pixel originated)
template <typename PixelIt>
class PixelMerger {
  using N = decltype(std::declval<PixelIt>()->count);
  struct Node {
    ThinPixel<N> pixel{};  // NOLINT
    std::size_t i{};       // NOLINT

    bool operator<(const Node &other) const noexcept;
    bool operator>(const Node &other) const noexcept;
    bool operator==(const Node &other) const noexcept;
    bool operator!=(const Node &other) const noexcept;
  };

  std::vector<PixelIt> _heads{};
  std::vector<PixelIt> _tails{};
  std::size_t _i{};

 public:
  using PixelT = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  class iterator;

  PixelMerger() = delete;
  PixelMerger(std::vector<PixelIt> head, std::vector<PixelIt> tail);

  template <typename ItOfPixelIt>
  PixelMerger(ItOfPixelIt first_head, ItOfPixelIt last_head, ItOfPixelIt first_tail);

  auto begin() const -> iterator;
  auto end() const noexcept -> iterator;

  [[nodiscard]] auto read_all() const -> std::vector<PixelT>;

  class iterator {
    PixelT _value{};

    using PQueueT = std::priority_queue<Node, std::vector<Node>, std::greater<>>;
    std::shared_ptr<PQueueT> _pqueue{};

    std::shared_ptr<std::vector<PixelIt>> _heads{};
    std::shared_ptr<std::vector<PixelIt>> _tails{};
    std::size_t _i{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = PixelT;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    explicit iterator(const std::vector<PixelIt> &heads, const std::vector<PixelIt> &tails);
    iterator(const iterator &other);
    iterator(iterator &&other) noexcept;
    ~iterator() noexcept = default;

    auto operator=(const iterator &other) -> iterator &;
    auto operator=(iterator &&other) noexcept -> iterator &;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    auto operator*() const noexcept -> const ThinPixel<N> &;
    auto operator->() const noexcept -> const ThinPixel<N> *;

    [[nodiscard]] auto operator++() -> iterator &;

   private:
    [[nodiscard]] auto next() -> PixelT;
    void replace_top_node();
  };
};
}  // namespace hictk::transformers

#include "./impl/pixel_merger_impl.hpp"  // NOLINT
