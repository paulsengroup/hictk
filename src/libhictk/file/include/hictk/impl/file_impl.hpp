// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/balancing/weights.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/pixel.hpp"

namespace hictk {

template <typename PixelSelectorT>
PixelSelector::PixelSelector(PixelSelectorT selector,
                             std::shared_ptr<const balancing::Weights> weights)
    : _sel(std::move(selector)), _weights(std::move(weights)) {
  assert(_weights);
}

template <typename N>
auto PixelSelector::begin([[maybe_unused]] bool sorted) const -> iterator<N> {
  return std::visit(
      [&](const auto& sel) {
        using T = std::decay_t<decltype(sel)>;
        if constexpr (std::is_same_v<cooler::PixelSelector, T>) {
          return iterator<N>{sel.template begin<N>(), sel.template end<N>()};
        } else {
          return iterator<N>{sel.template begin<N>(sorted), sel.template end<N>()};
        }
      },
      _sel);
}

template <typename N>
auto PixelSelector::end() const -> iterator<N> {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) { return iterator<N>{sel.template end<N>(), sel.template end<N>()}; },
      _sel);
}

template <typename N>
auto PixelSelector::cbegin(bool sorted) const -> iterator<N> {
  return begin<N>(sorted);
}

template <typename N>
auto PixelSelector::cend() const -> iterator<N> {
  return end<N>();
}

template <typename N>
std::vector<Pixel<N>> PixelSelector::read_all() const {
  assert(!_sel.valueless_by_exception());
  return std::visit([&](const auto& sel) { return sel.template read_all<N>(); }, _sel);
}

template <typename PixelSelectorT>
constexpr const PixelSelectorT& PixelSelector::get() const {
  assert(!_sel.valueless_by_exception());
  return std::get<PixelSelectorT>(_sel);
}

template <typename PixelSelectorT>
constexpr PixelSelectorT& PixelSelector::get() {
  assert(!_sel.valueless_by_exception());
  return std::get<PixelSelectorT>(_sel);
}

constexpr auto PixelSelector::get() const noexcept -> const PixelSelectorVar& { return _sel; }
constexpr auto PixelSelector::get() noexcept -> PixelSelectorVar& { return _sel; }

template <typename N>
template <typename It>
PixelSelector::iterator<N>::iterator(It it, It end)
    : _it(std::move(it)), _sentinel(std::move(end)) {}

template <typename N>  // NOLINTNEXTLINE(bugprone-exception-escape)
bool PixelSelector::iterator<N>::operator==(const iterator& other) const noexcept {
  assert(!_it.valueless_by_exception());
  return operator_eq(_it, other._it);
}

template <typename N>
bool PixelSelector::iterator<N>::operator!=(const iterator& other) const noexcept {
  return !(*this == other);
}

template <typename N>
auto PixelSelector::iterator<N>::operator*() const -> const_reference {
  assert(!_it.valueless_by_exception());
  return std::visit([&](const auto& it) -> const_reference { return *it; }, _it);
}

template <typename N>
auto PixelSelector::iterator<N>::operator->() const -> const_pointer {
  assert(!_it.valueless_by_exception());
  return std::visit([&](const auto& it) -> const_pointer { return &*it; }, _it);
}

template <typename N>
auto PixelSelector::iterator<N>::operator++() -> iterator& {
  assert(!_it.valueless_by_exception());
  std::visit([&](auto& it) { ++it; }, _it);
  return *this;
}

template <typename N>
auto PixelSelector::iterator<N>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N>
template <typename IteratorT>
[[nodiscard]] constexpr const IteratorT& PixelSelector::iterator<N>::get() const {
  return std::get<IteratorT>(_it);
}

template <typename N>
template <typename IteratorT>
[[nodiscard]] constexpr IteratorT& PixelSelector::iterator<N>::get() {
  return std::get<IteratorT>(_it);
}

template <typename N>
constexpr auto PixelSelector::iterator<N>::get() const noexcept -> const IteratorVar& {
  return _it;
}

template <typename N>
constexpr auto PixelSelector::iterator<N>::get() noexcept -> IteratorVar& {
  return _it;
}

template <typename N>  // NOLINTNEXTLINE(bugprone-exception-escape)
bool PixelSelector::iterator<N>::operator_eq(const IteratorVar& itv1,
                                             const IteratorVar& itv2) noexcept {
  assert(!itv1.valueless_by_exception());
  assert(!itv2.valueless_by_exception());
  return std::visit(
      [&](const auto& it1) {
        using T = std::decay_t<decltype(it1)>;
        const auto* it2 = std::get_if<T>(&itv2);
        return !!it2 && it1 == *it2;
      },
      itv1);
}

template <typename N>  // NOLINTNEXTLINE(bugprone-exception-escape)
bool PixelSelector::iterator<N>::operator_neq(const IteratorVar& itv1,
                                              const IteratorVar& itv2) noexcept {
  return !(itv1 == itv2);
}

constexpr bool File::is_hic() const noexcept {
  assert(!_fp.valueless_by_exception());
  return std::holds_alternative<hic::File>(_fp);
}

constexpr bool File::is_cooler() const noexcept { return !is_hic(); }

template <typename FileT>
constexpr const FileT& File::get() const {
  return std::get<FileT>(_fp);
}

template <typename FileT>
constexpr FileT& File::get() {
  return std::get<FileT>(_fp);
}

constexpr auto File::get() const noexcept -> const FileVar& { return _fp; }
constexpr auto File::get() noexcept -> FileVar& { return _fp; }

}  // namespace hictk
