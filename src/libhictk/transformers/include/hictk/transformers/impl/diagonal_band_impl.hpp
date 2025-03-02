// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

// NOLINTNEXTLINE(*-concat-nested-namespaces)
namespace hictk::transformers {

namespace internal {
template <typename T, typename = std::void_t<>>
inline constexpr bool has_bin1_id_member = false;

template <typename T>
inline constexpr bool has_bin1_id_member<T, std::void_t<decltype(std::declval<T>().bin1_id)>> =
    true;

template <typename T, typename = std::void_t<>>
inline constexpr bool has_jump_to_next_row_member_fx = false;

template <typename T>
inline constexpr bool
    has_jump_to_next_row_member_fx<T, std::void_t<decltype(std::declval<T>().jump_to_next_row())>> =
        true;
}  // namespace internal

template <typename PixelIt>
inline DiagonalBand<PixelIt>::DiagonalBand(PixelIt first, PixelIt last,
                                           std::uint64_t num_bins) noexcept
    : _first(std::move(first)), _last(std::move(last)), _num_bins(num_bins) {
  if (num_bins == 0) {
    _first = _last;
  }
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::begin() const -> iterator {
  return iterator(_first, _last, _num_bins);
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::end() const -> iterator {
  return iterator::at_end(_last);
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::cbegin() const -> iterator {
  return begin();
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::cend() const -> iterator {
  return end();
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::read_all() const -> std::vector<Pixel> {
  // We push_back into buff to avoid traversing pixels twice (once to figure out the vector size,
  // and a second time to copy the actual data)
  std::vector<Pixel> buff{};
  std::copy(begin(), end(), std::back_inserter(buff));
  return buff;
}

template <typename PixelIt>
inline DiagonalBand<PixelIt>::iterator::iterator(PixelIt first, PixelIt last,
                                                 std::uint64_t num_bins) noexcept
    : _it(std::move(first)), _last(std::move(last)), _num_bins(num_bins) {}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::iterator::at_end(const PixelIt &it) noexcept -> iterator {
  return {it, it, 0};
}

template <typename PixelIt>
inline bool DiagonalBand<PixelIt>::iterator::operator==(const iterator &other) const noexcept {
  return _it == other._it;
}

template <typename PixelIt>
inline bool DiagonalBand<PixelIt>::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::iterator::operator*() const -> const_reference {
  assert(_it != _last);
  return *_it;
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::iterator::operator->() const -> const_pointer {
  return _it.operator->();
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::iterator::operator++() -> iterator & {
  if (++_it == _last) {
    return *this;
  }

  const auto p = *_it;

  if constexpr (internal::has_bin1_id_member<Pixel>) {
    if (p.bin2_id - p.bin1_id < _num_bins) {
      return *this;
    }
  } else {
    if (p.bin2.id() - p.bin1.id() < _num_bins) {
      return *this;
    }
  }

  if constexpr (internal::has_jump_to_next_row_member_fx<PixelIt>) {
    _it->jump_to_next_row();
    return *this;
  }
  return ++(*this);
}

template <typename PixelIt>
inline auto DiagonalBand<PixelIt>::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

}  // namespace hictk::transformers
