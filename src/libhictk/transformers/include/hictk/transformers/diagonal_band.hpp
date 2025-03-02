// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <type_traits>
#include <vector>

#include "hictk/type_traits.hpp"

namespace hictk::transformers {

template <typename PixelIt>
class DiagonalBand {
  using Pixel = remove_cvref_t<decltype(*std::declval<PixelIt>())>;
  PixelIt _first{};
  PixelIt _last{};
  std::uint64_t _num_bins{};

 public:
  class iterator;

  DiagonalBand(PixelIt first, PixelIt last, std::uint64_t num_bins) noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;

  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  [[nodiscard]] auto read_all() const -> std::vector<Pixel>;

  class iterator {
    PixelIt _it{};
    PixelIt _last{};
    std::uint64_t _num_bins{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Pixel;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(PixelIt first, PixelIt last, std::uint64_t num_bins) noexcept;
    static auto at_end(const PixelIt &it) noexcept -> iterator;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
  };
};

}  // namespace hictk::transformers

#include "./impl/diagonal_band_impl.hpp"
