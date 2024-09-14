// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::balancing {

inline Weights::Weights(std::vector<double> weights, Type type)
    : _weights(std::make_shared<WeightVect>(std::move(weights))), _type(type) {
  if (_type != Type::MULTIPLICATIVE && _type != Type::DIVISIVE) {
    throw std::runtime_error("Weight type must be either MULTIPLICATIVE or DIVISIVE");
  }
}

inline Weights::Weights(std::vector<double> weights, std::string_view name)
    : _weights(std::make_shared<WeightVect>(std::move(weights))), _type(Weights::infer_type(name)) {
  assert(_type != Type::INFER);
  if (_type == Type::UNKNOWN) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), name));
  }
}

inline Weights::Weights(double weight, std::size_t size, Type type)
    : _weights(ConstWeight{weight, size}), _type(type) {
  if (_type != Type::MULTIPLICATIVE && _type != Type::DIVISIVE) {
    throw std::runtime_error("Weight type must be either MULTIPLICATIVE or DIVISIVE");
  }
}

inline Weights::Weights(double weight, std::size_t size, std::string_view name)
    : _weights(ConstWeight{weight, size}), _type(Weights::infer_type(name)) {
  assert(_type != Type::INFER);
  if (_type == Type::UNKNOWN) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), name));
  }
}

inline Weights::Weights(std::variant<ConstWeight, WeightVectPtr> weights, Type type_) noexcept
    : _weights(std::move(weights)), _type(type_) {
  assert(_type == Type::MULTIPLICATIVE || _type == Type::DIVISIVE);
}

inline Weights::operator bool() const noexcept { return !empty(); }

inline double Weights::operator[](std::size_t i) const noexcept {
  if (is_constant()) {
    return std::get<ConstWeight>(_weights).w;
  }

  assert(std::get<WeightVectPtr>(_weights));
  return (*std::get<WeightVectPtr>(_weights))[i];
}

inline double Weights::at(std::size_t i) const {
  if (is_constant()) {
    const auto& [w, size] = std::get<ConstWeight>(_weights);

    if (i >= size) {
      throw std::out_of_range("Weights::at()");
    }

    return w;
  }

  return std::get<WeightVectPtr>(_weights)->at(i);
}

inline double Weights::at(std::size_t i, Type type_) const {
  if (HICTK_UNLIKELY(type_ != Type::MULTIPLICATIVE && type_ != Type::DIVISIVE)) {
    throw std::logic_error("Type should be Type::MULTIPLICATIVE or Type::DIVISIVE");
  }

  if (type_ == type()) {
    return at(i);
  }

  return 1.0 / at(i);
}

inline auto Weights::begin(Type type_) const -> iterator {
  if (type_ == Type::UNKNOWN) {
    throw std::logic_error("Weights::begin(): type cannot be UNKNOWN");
  }
  return cbegin(type_);
}

inline auto Weights::end(Type type_) const -> iterator {
  if (type_ == Type::UNKNOWN) {
    throw std::logic_error("Weights::end(): type cannot be UNKNOWN");
  }
  return cend(type_);
}

inline auto Weights::cbegin(Type type_) const -> iterator {
  if (type_ == Type::UNKNOWN) {
    throw std::logic_error("Weights::cbegin(): type cannot be UNKNOWN");
  }
  if (type_ == Type::INFER) {
    type_ = type();
  }

  const auto reciprocal = type_ != type();

  if (is_constant()) {
    return {std::get<ConstWeight>(_weights), 0, reciprocal};
  }

  assert(std::get<WeightVectPtr>(_weights));
  return {std::get<WeightVectPtr>(_weights)->begin(), reciprocal};
}

inline auto Weights::cend(Type type_) const -> iterator {
  if (type_ == Type::UNKNOWN) {
    throw std::logic_error("Weights::cend(): type cannot be UNKNOWN");
  }

  if (type_ == Type::INFER) {
    type_ = type();
  }
  const auto reciprocal = type_ != type();

  if (is_constant()) {
    const auto& weights = std::get<ConstWeight>(_weights);
    return iterator{weights, weights.size, reciprocal};
  }

  assert(std::get<WeightVectPtr>(_weights));
  return iterator{std::get<WeightVectPtr>(_weights)->end(), reciprocal};
}

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

inline Weights Weights::operator()(Type type_) const {
  if (HICTK_UNLIKELY(type_ != Type::MULTIPLICATIVE && type_ != Type::DIVISIVE)) {
    throw std::logic_error("Type should be Type::MULTIPLICATIVE or Type::DIVISIVE");
  }

  if (type_ == type()) {
    return *this;
  }

  if (is_constant()) {
    const auto& [w, size] = std::get<ConstWeight>(_weights);
    return {ConstWeight{1.0 / w, size}, type_};
  }

  const auto& buff = std::get<WeightVectPtr>(_weights);
  assert(!!buff);

  auto weights = std::make_shared<WeightVect>(buff->size());
  std::transform(buff->begin(), buff->end(), weights->begin(),
                 [](const auto n) { return 1.0 / n; });

  return {std::move(weights), type_};
}

constexpr auto Weights::type() const noexcept -> Type { return _type; }

inline std::size_t Weights::size() const noexcept {
  if (is_constant()) {
    return std::get<ConstWeight>(_weights).size;
  }
  const auto& weights = std::get<WeightVectPtr>(_weights);
  if (!weights) {
    return 0;
  }
  return weights->size();
}

inline bool Weights::empty() const noexcept { return size() == 0; }

inline auto Weights::infer_type(std::string_view name) noexcept -> Type {
  constexpr std::array<std::pair<std::string_view, Type>, 14> mappings{
      {{"VC", Type::DIVISIVE},
       {"INTER_VC", Type::DIVISIVE},
       {"GW_VC", Type::DIVISIVE},
       {"VC_SQRT", Type::DIVISIVE},
       {"KR", Type::DIVISIVE},
       {"INTER_KR", Type::DIVISIVE},
       {"GW_KR", Type::DIVISIVE},
       {"SCALE", Type::DIVISIVE},
       {"INTER_SCALE", Type::DIVISIVE},
       {"GW_SCALE", Type::DIVISIVE},
       {"ICE", Type::MULTIPLICATIVE},
       {"INTER_ICE", Type::MULTIPLICATIVE},
       {"GW_ICE", Type::MULTIPLICATIVE},
       {"weight", Type::MULTIPLICATIVE}}};

  auto it = std::find_if(mappings.begin(), mappings.end(),
                         [&](const auto& p) { return p.first == name; });
  if (it == mappings.end()) {
    return Weights::Type::UNKNOWN;
  }
  return it->second;
}

inline void Weights::rescale(double scaling_factor) noexcept {
  if (is_constant()) {
    auto& w = std::get<ConstWeight>(_weights).w;
    w *= std::sqrt(scaling_factor);
  } else {
    auto& weights = std::get<WeightVectPtr>(_weights);
    assert(!!weights);
    if (weights.use_count() != 1) {
      _weights = std::make_shared<WeightVect>(*weights);
    }
    std::transform(weights->begin(), weights->end(), weights->begin(),
                   [&](auto w) { return w * std::sqrt(scaling_factor); });
  }
}

inline void Weights::rescale(const std::vector<double>& scaling_factors,
                             const std::vector<std::uint64_t>& offsets) {
  if (scaling_factors.empty()) {
    throw std::runtime_error("scaling_factors cannot be empty");
  }

  if (scaling_factors.size() + 1 != offsets.size()) {
    throw std::runtime_error("offsets must have a size of scaling_factors.size() + 1");
  }

  if (offsets.front() != 0) {
    throw std::runtime_error("the first offset should be 0");
  }

  if (offsets.back() != size()) {
    throw std::runtime_error("the last offset should be the size of the weight vector");
  }

  if (is_constant()) {
    if (scaling_factors.size() == 1) {
      rescale(scaling_factors.front());
      return;
    }
    throw std::runtime_error(
        "rescaling ConstWeight with multiple scaling factors is not supported");
  }

  if (scaling_factors.size() == 1) {
    rescale(scaling_factors.front());
    return;
  }

  if (!std::is_sorted(offsets.begin(), offsets.end())) {
    throw std::runtime_error("offset vector is not sorted in ascending order");
  }

  const auto& weights = std::get<WeightVectPtr>(_weights);
  assert(!!weights);

  if (weights.use_count() != 1) {
    _weights = std::make_shared<WeightVect>(*weights);
  }

  for (std::size_t i = 0; i < scaling_factors.size(); ++i) {
    auto first = weights->begin() + std::ptrdiff_t(offsets[i]);
    auto last = weights->begin() + std::ptrdiff_t(offsets[i + 1]);
    std::transform(first, last, first,
                   [s = scaling_factors[i]](const double w) { return w * std::sqrt(s); });
  }
}

inline std::vector<double> Weights::to_vector(Type type_) const {
  assert(type_ != Type::UNKNOWN);
  if (type_ == Type::INFER) {
    type_ = type();
  }

  const auto weights = (*this)(type_);

  if (is_constant()) {
    const auto& [w, size] = std::get<ConstWeight>(weights._weights);
    return std::vector<double>(size, w);
  }

  assert(std::get<WeightVectPtr>(weights._weights));
  return *std::get<WeightVectPtr>(weights._weights);
}

inline bool Weights::is_constant() const noexcept {
  return std::holds_alternative<ConstWeight>(_weights);
}

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

inline auto Weights::iterator::ConstIt::operator*() const noexcept(ndebug_defined()) -> double {
  assert(!!value);
  bound_check(0, false);
  return value->w;
}

inline auto Weights::iterator::ConstIt::operator[]([[maybe_unused]] std::ptrdiff_t i_) const
    noexcept(ndebug_defined()) -> double {
  bound_check(i_, false);
  return **this;
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

constexpr void Weights::iterator::ConstIt::bound_check(std::ptrdiff_t offset, bool end_ok) const {
  if constexpr (ndebug_not_defined()) {
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
        throw std::logic_error(fmt::format(FMT_STRING("Invalid offset {}: {} + {} >= {}"), offset,
                                           i, offset, value->size));
      }
    }
  }
}

inline Weights::iterator::iterator(std::vector<double>::const_iterator it, bool reciprocal)
    : _it(std::move(it)), _reciprocal(reciprocal) {}

inline Weights::iterator::iterator(const ConstWeight& weight, std::size_t offset, bool reciprocal)
    : _it(ConstIt{&weight, static_cast<std::ptrdiff_t>(offset)}), _reciprocal(reciprocal) {
  assert(offset <= weight.size);
}

inline Weights::iterator::iterator(ConstIt it, bool reciprocal) noexcept
    : _it(std::move(it)), _reciprocal(reciprocal) {}

inline bool Weights::iterator::is_constant() const noexcept {
  return std::holds_alternative<ConstIt>(_it);
}

inline bool Weights::iterator::operator==(const iterator& other) const noexcept {
  if (_reciprocal != other._reciprocal) {
    return false;
  }

  if (HICTK_LIKELY(is_constant() == other.is_constant())) {
    return std::visit(
        [&](const auto& it1) {
          using T = remove_cvref_t<decltype(it1)>;
          const auto& it2 = std::get<T>(other._it);
          return it1 == it2;
        },
        _it);
  }
  return false;
}

inline bool Weights::iterator::operator!=(const iterator& other) const noexcept {
  return !(*this == other);
}

inline bool Weights::iterator::operator<(const iterator& other) const {
  assert(_reciprocal == other._reciprocal);

  if (HICTK_LIKELY(is_constant() == other.is_constant())) {
    return std::visit(
        [&](const auto& it1) {
          using T = remove_cvref_t<decltype(it1)>;
          const auto& it2 = std::get<T>(other._it);
          return it1 < it2;
        },
        _it);
  }
  throw std::logic_error("caught attempt to compare iterators of different type");
}

inline bool Weights::iterator::operator<=(const iterator& other) const {
  assert(_reciprocal == other._reciprocal);

  if (HICTK_LIKELY(is_constant() == other.is_constant())) {
    return std::visit(
        [&](const auto& it1) {
          using T = remove_cvref_t<decltype(it1)>;
          const auto& it2 = std::get<T>(other._it);
          return it1 <= it2;
        },
        _it);
  }
  throw std::logic_error("caught attempt to compare iterators of different type");
}

inline bool Weights::iterator::operator>(const iterator& other) const {
  assert(_reciprocal == other._reciprocal);

  if (HICTK_LIKELY(is_constant() == other.is_constant())) {
    return std::visit(
        [&](const auto& it1) {
          using T = remove_cvref_t<decltype(it1)>;
          const auto& it2 = std::get<T>(other._it);
          return it1 > it2;
        },
        _it);
  }
  throw std::logic_error("caught attempt to compare iterators of different type");
}

inline bool Weights::iterator::operator>=(const iterator& other) const {
  assert(_reciprocal == other._reciprocal);

  if (HICTK_LIKELY(is_constant() == other.is_constant())) {
    return std::visit(
        [&](const auto& it1) {
          using T = remove_cvref_t<decltype(it1)>;
          const auto& it2 = std::get<T>(other._it);
          return it1 >= it2;
        },
        _it);
  }
  throw std::logic_error("caught attempt to compare iterators of different type");
}

inline auto Weights::iterator::operator*() const -> value_type {
  const auto w = std::visit([&](const auto& it) { return *it; }, _it);
  if (_reciprocal) {
    return 1.0 / w;
  }

  return w;
}

inline auto Weights::iterator::operator[](difference_type i) const -> value_type {
  const auto w = std::visit([&](const auto& it) { return it[i]; }, _it);
  if (_reciprocal) {
    return 1.0 / w;
  }

  return w;
}

inline auto Weights::iterator::operator++() -> iterator& {
  std::visit([&](auto& it) { ++it; }, _it);
  return *this;
}

inline auto Weights::iterator::operator++(int) -> iterator {
  return std::visit([&](auto& it_) -> iterator { return {it_++, _reciprocal}; }, _it);
}

inline auto Weights::iterator::operator+=(difference_type i) -> iterator& {
  std::visit([&](auto& it) { it += i; }, _it);
  return *this;
}

inline auto Weights::iterator::operator+(difference_type i) const -> iterator {
  return std::visit([&](const auto& it_) -> iterator { return {it_ + i, _reciprocal}; }, _it);
}

inline auto Weights::iterator::operator--() -> iterator& {
  std::visit([&](auto& it) { --it; }, _it);
  return *this;
}

inline auto Weights::iterator::operator--(int) -> iterator {
  return std::visit([&](auto& it_) -> iterator { return {it_--, _reciprocal}; }, _it);
}

inline auto Weights::iterator::operator-=(difference_type i) -> iterator& {
  std::visit([&](auto& it) { it -= i; }, _it);
  return *this;
}

inline auto Weights::iterator::operator-(difference_type i) const -> iterator {
  return std::visit([&](const auto& it_) -> iterator { return {it_ - i, _reciprocal}; }, _it);
}

inline auto Weights::iterator::operator-(const iterator& other) const -> difference_type {
  assert(_reciprocal == other._reciprocal);

  if (HICTK_LIKELY(is_constant() == other.is_constant())) {
    return std::visit(
        [&](const auto& it1) {
          using T = remove_cvref_t<decltype(it1)>;
          const auto& it2 = std::get<T>(other._it);
          return it1 - it2;
        },
        _it);
  }
  throw std::logic_error("caught attempt to compare iterators of different type");
}

}  // namespace hictk::balancing
