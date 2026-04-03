// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifndef NDEBUG
#include <fmt/format.h>
#endif

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::balancing {

template <typename N>
inline ThinPixel<N> Weights::balance(ThinPixel<N> p) const {
  p.count = balance<N>(p.bin1_id, p.bin2_id, p.count);
  return p;
}

template <typename N>
inline Pixel<N> Weights::balance(Pixel<N> p) const {
  p.count = balance<N>(p.coords.bin1().id(), p.coords.bin2().id(), p.count);
  return p;
}

template <typename N1, typename N2>
inline N1 Weights::balance(std::uint64_t bin1_id, std::uint64_t bin2_id, N2 count) const {
  assert(std::is_floating_point_v<N1>);
  const auto w1 = at(conditional_static_cast<std::size_t>(bin1_id));
  const auto w2 = at(conditional_static_cast<std::size_t>(bin2_id));

  auto count_ = conditional_static_cast<double>(count);

  if (type() == Weights::Type::MULTIPLICATIVE) {
    count_ *= w1 * w2;
  } else {
    assert(type() == Weights::Type::DIVISIVE);
    count_ /= w1 * w2;
  }
  return conditional_static_cast<N1>(count_);
}

constexpr auto Weights::type() const noexcept -> Type { return _type; }

constexpr bool Weights::iterator::ConstIt::operator==(const ConstIt& other) const noexcept {
  return value == other.value && i == other.i;
}

constexpr bool Weights::iterator::ConstIt::operator!=(const ConstIt& other) const noexcept {
  return !(*this == other);
}

constexpr bool Weights::iterator::ConstIt::operator<(const ConstIt& other) const noexcept {
  assert(value == other.value);
  return i < other.i;
}

constexpr bool Weights::iterator::ConstIt::operator<=(const ConstIt& other) const noexcept {
  assert(value == other.value);
  return i <= other.i;
}

constexpr bool Weights::iterator::ConstIt::operator>(const ConstIt& other) const noexcept {
  assert(value == other.value);
  return i > other.i;
}

constexpr bool Weights::iterator::ConstIt::operator>=(const ConstIt& other) const noexcept {
  assert(value == other.value);
  return i >= other.i;
}

constexpr auto Weights::iterator::ConstIt::operator++() noexcept(ndebug_defined()) -> ConstIt& {
  bound_check(1);
  ++i;
  return *this;
}

constexpr auto Weights::iterator::ConstIt::operator++(int) noexcept(ndebug_defined()) -> ConstIt {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

constexpr auto Weights::iterator::ConstIt::operator+=(std::ptrdiff_t i_) noexcept(ndebug_defined())
    -> ConstIt& {
  bound_check(i_);
  i += i_;
  return *this;
}

constexpr auto Weights::iterator::ConstIt::operator+(std::ptrdiff_t i_) const
    noexcept(ndebug_defined()) -> ConstIt {
  auto it = *this;
  return it += i_;
}

constexpr auto Weights::iterator::ConstIt::operator--() noexcept(ndebug_defined()) -> ConstIt& {
  bound_check(-1);
  --i;
  return *this;
}

constexpr auto Weights::iterator::ConstIt::operator--(int) noexcept(ndebug_defined()) -> ConstIt {
  bound_check(-1);
  auto it = *this;
  std::ignore = --(*this);
  return it;
}

constexpr auto Weights::iterator::ConstIt::operator-=(std::ptrdiff_t i_) noexcept(ndebug_defined())
    -> ConstIt& {
  bound_check(i_);
  i -= i_;
  return *this;
}

constexpr auto Weights::iterator::ConstIt::operator-(std::ptrdiff_t i_) const
    noexcept(ndebug_defined()) -> ConstIt {
  auto it = *this;
  return it -= i_;
}

constexpr auto Weights::iterator::ConstIt::operator-(const ConstIt& other) const
    noexcept(ndebug_defined()) -> std::ptrdiff_t {
  assert(value == other.value);
  return i - other.i;
}

constexpr void Weights::iterator::ConstIt::bound_check([[maybe_unused]] std::ptrdiff_t offset,
                                                       [[maybe_unused]] bool end_ok) const {
#ifndef NDEBUG  // GCC8 does not like it when we use if constexpr in this context
  assert(!!value);
  if (offset < 0 && offset > i) {
    throw std::logic_error(
        fmt::format(FMT_STRING("Invalid offset {}: {} - {} < 0"), offset, i, offset));
  }

  if (end_ok) {
    if (i + offset > static_cast<std::ptrdiff_t>(value->size)) {
      throw std::logic_error(fmt::format(FMT_STRING("Invalid offset {}: {} + {} > {}"), offset, i,
                                         offset, value->size));
    }
  } else {
    if (i + offset > static_cast<std::ptrdiff_t>(value->size)) {
      throw std::logic_error(fmt::format(FMT_STRING("Invalid offset {}: {} + {} >= {}"), offset, i,
                                         offset, value->size));
    }
  }
#endif
}

}  // namespace hictk::balancing
