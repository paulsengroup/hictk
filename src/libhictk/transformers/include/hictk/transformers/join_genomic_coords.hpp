// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <iterator>
#include <memory>
#include <type_traits>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename PixelIt>
class JoinGenomicCoords {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using N = decltype(std::declval<PixelT>().count);
  PixelIt _first{};
  PixelIt _last{};
  std::shared_ptr<const BinTable> _bins{};

 public:
  class iterator;

  JoinGenomicCoords(PixelIt first, PixelIt last, std::shared_ptr<const BinTable> bins);

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;

  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  [[nodiscard]] auto read_all() const -> std::vector<Pixel<N>>;

  class iterator {
    PixelIt _it{};
    std::shared_ptr<const BinTable> _bins{};
    mutable Pixel<N> _value{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Pixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(PixelIt it, std::shared_ptr<const BinTable> bins);
    static auto at_end(PixelIt it, std::shared_ptr<const BinTable> bins) -> iterator;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
  };
};

}  // namespace hictk::transformers

#include "./impl/join_genomic_coords_impl.hpp"  // NOLINT
