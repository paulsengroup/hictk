// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>
#include <variant>

#include "hictk/balancing/methods.hpp"

namespace hictk::cooler {

template <typename T>
inline bool File::has_pixel_of_type() const noexcept {
  return std::holds_alternative<T>(_pixel_variant);
}

template <typename N>
inline PixelSelector::iterator<N> File::begin(std::string_view weight_name) const {
  return fetch(normalization_ptr(balancing::Method(weight_name))).template begin<N>();
}

template <typename N>
inline PixelSelector::iterator<N> File::cbegin(std::string_view weight_name) const {
  return begin<N>(weight_name);
}

template <typename N>
inline PixelSelector::iterator<N> File::end(std::string_view weight_name) const {
  return fetch(normalization_ptr(balancing::Method(weight_name))).template end<N>();
}

template <typename N>
inline PixelSelector::iterator<N> File::cend(std::string_view weight_name) const {
  return end<N>(weight_name);
}

}  // namespace hictk::cooler
