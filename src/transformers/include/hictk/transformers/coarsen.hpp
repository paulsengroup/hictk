// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>

#include <cstdint>
#include <memory>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename PixelIt>
class CoarsenPixels {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using N = decltype(std::declval<PixelT>().count);

  PixelIt _first{};
  PixelIt _last{};
  std::size_t _factor{};

 public:
  class iterator;

  CoarsenPixels(PixelIt first_pixel, PixelIt last_pixel, std::size_t factor);

  auto begin() const -> iterator;
  auto end() const -> iterator;

  auto cbegin() const -> iterator;
  auto cend() const -> iterator;

  [[nodiscard]] auto read_all() const -> std::vector<ThinPixel<N>>;

  class iterator {
    struct PixelCmp {
      [[nodiscard]] inline bool operator()(const ThinPixel<N> &p1,
                                           const ThinPixel<N> &p2) const noexcept;
    };

    using PixelMerger = phmap::btree_set<ThinPixel<N>, PixelCmp>;
    using It = typename PixelMerger::const_iterator;

    PixelIt _pixel_it{};
    PixelIt _pixel_last{};
    std::shared_ptr<PixelMerger> _merger{};
    It _it{};

    std::uint64_t _bin1_id_end{};
    std::size_t _factor{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = ThinPixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(PixelIt first, PixelIt last, std::size_t factor);
    static auto at_end(PixelIt last, std::size_t factor) -> iterator;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    // Current implementation does not allow for an efficient implementation of it++
    auto operator++(int) -> iterator;

   private:
    void process_next_row();
  };
};

}  // namespace hictk::transformers

#include "../../../coarsen_impl.hpp"
