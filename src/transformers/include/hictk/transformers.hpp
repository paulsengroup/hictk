// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/hic/pixel_selector.hpp"

namespace hictk::transformers {

class JoinGenomicCoords {
  std::variant<cooler::PixelSelector<>, hic::PixelSelector> _sel;
  std::shared_ptr<const BinTable> _bins{};

 public:
  template <typename N>
  class iterator;

  JoinGenomicCoords(cooler::PixelSelector<> sel, std::shared_ptr<const BinTable> bins);
  JoinGenomicCoords(hic::PixelSelector sel, std::shared_ptr<const BinTable> bins);

  template <typename N>
  auto begin() -> iterator<N>;
  template <typename N>
  auto end() -> iterator<N>;

  template <typename N>
  auto cbegin() -> iterator<N>;
  template <typename N>
  auto cend() -> iterator<N>;

  template <typename N>
  [[nodiscard]] std::vector<Pixel<N>> read_all() const;

  template <typename N>
  class iterator {
    using CoolerIt = cooler::PixelSelector<>::iterator<N>;
    using HicIt = hic::PixelSelector::iterator<N>;
    std::variant<CoolerIt, HicIt> _it{CoolerIt{}};
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

    explicit iterator(CoolerIt it, std::shared_ptr<const BinTable> bins);
    explicit iterator(HicIt it, std::shared_ptr<const BinTable> bins);

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
  };
};

}  // namespace hictk::transformers

#include "../../transformers_impl.hpp"
