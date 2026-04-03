// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <optional>

namespace hictk::cooler {

template <typename T>
inline auto Dataset::cbegin(std::optional<std::ptrdiff_t> chunk_size) const -> iterator<T> {
  return iterator<T>(*this, chunk_size);
}

template <typename T>
inline auto Dataset::cend(std::optional<std::ptrdiff_t> chunk_size) const -> iterator<T> {
  return iterator<T>::make_end_iterator(*this, chunk_size);
}

template <typename T>
inline auto Dataset::begin(std::optional<std::ptrdiff_t> chunk_size) const -> iterator<T> {
  return cbegin<T>(chunk_size);
}

template <typename T>
inline auto Dataset::end(std::optional<std::ptrdiff_t> chunk_size) const -> iterator<T> {
  return cend<T>(chunk_size);
}

}  // namespace hictk::cooler
