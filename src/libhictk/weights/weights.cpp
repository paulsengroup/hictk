// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/weights.hpp"

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
#include <utility>
#include <variant>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk {

Weights::Weights(std::vector<double> weights, Type type)
    : _weights(std::make_shared<WeightVect>(std::move(weights))), _type(type) {
  if (_type != Type::MULTIPLICATIVE && _type != Type::DIVISIVE) {
    throw std::runtime_error("Weight type must be either MULTIPLICATIVE or DIVISIVE");
  }
}

Weights::Weights(std::vector<double> weights, std::string_view name)
    : _weights(std::make_shared<WeightVect>(std::move(weights))), _type(Weights::infer_type(name)) {
  assert(_type != Type::INFER);
  if (_type == Type::UNKNOWN) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), name));
  }
}

Weights::Weights(double weight, std::size_t size, Type type)
    : _weights(ConstWeight{weight, size}), _type(type) {
  if (_type != Type::MULTIPLICATIVE && _type != Type::DIVISIVE) {
    throw std::runtime_error("Weight type must be either MULTIPLICATIVE or DIVISIVE");
  }
}

Weights::Weights(double weight, std::size_t size, std::string_view name)
    : _weights(ConstWeight{weight, size}), _type(infer_type(name)) {
  assert(_type != Type::INFER);
  if (_type == Type::UNKNOWN) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), name));
  }
}

Weights::Weights(std::variant<ConstWeight, WeightVectPtr> weights, Type type_) noexcept
    : _weights(std::move(weights)), _type(type_) {
  assert(_type == Type::MULTIPLICATIVE || _type == Type::DIVISIVE);
}

Weights::operator bool() const noexcept { return !empty(); }

// NOLINTNEXTLINE(bugprone-exception-escape)
double Weights::operator[](std::size_t i) const noexcept {
  assert(!_weights.valueless_by_exception());
  if (is_constant()) {
    return std::get<ConstWeight>(_weights).w;
  }

  assert(std::get<WeightVectPtr>(_weights));
  return (*std::get<WeightVectPtr>(_weights))[i];
}

double Weights::at(std::size_t i) const {
  assert(!_weights.valueless_by_exception());
  if (is_constant()) {
    const auto& [w, size] = std::get<ConstWeight>(_weights);

    if (i >= size) {
      throw std::out_of_range("Weights::at()");
    }

    return w;
  }

  return std::get<WeightVectPtr>(_weights)->at(i);
}

double Weights::at(std::size_t i, Type type_) const {
  if (HICTK_UNLIKELY(type_ != Type::MULTIPLICATIVE && type_ != Type::DIVISIVE)) {
    throw std::logic_error("Type should be Type::MULTIPLICATIVE or Type::DIVISIVE");
  }

  if (type_ == type()) {
    return at(i);
  }

  return 1.0 / at(i);
}

auto Weights::begin(Type type_) const -> iterator {
  if (type_ == Type::UNKNOWN) {
    throw std::logic_error("Weights::begin(): type cannot be UNKNOWN");
  }
  return cbegin(type_);
}

auto Weights::end(Type type_) const -> iterator {
  if (type_ == Type::UNKNOWN) {
    throw std::logic_error("Weights::end(): type cannot be UNKNOWN");
  }
  return cend(type_);
}

auto Weights::cbegin(Type type_) const -> iterator {
  assert(!_weights.valueless_by_exception());
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

auto Weights::cend(Type type_) const -> iterator {
  assert(!_weights.valueless_by_exception());
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

Weights Weights::operator()(Type type_) const {
  assert(!_weights.valueless_by_exception());
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

// NOLINTNEXTLINE(bugprone-exception-escape)
std::size_t Weights::size() const noexcept {
  assert(!_weights.valueless_by_exception());
  if (is_constant()) {
    return std::get<ConstWeight>(_weights).size;
  }
  const auto& weights = std::get<WeightVectPtr>(_weights);
  if (!weights) {
    return 0;
  }
  return weights->size();
}

bool Weights::empty() const noexcept { return size() == 0; }

auto Weights::infer_type(std::string_view name) noexcept -> Type {
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

  // NOLINTNEXTLINE(*-qualified-auto)
  const auto it = std::find_if(mappings.begin(), mappings.end(),
                               [&](const auto& p) { return p.first == name; });
  if (it == mappings.end()) {
    return Weights::Type::UNKNOWN;
  }
  return it->second;
}

void Weights::rescale(double scaling_factor) {
  assert(!_weights.valueless_by_exception());
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

void Weights::rescale(const std::vector<double>& scaling_factors,
                      const std::vector<std::uint64_t>& offsets) {
  assert(!_weights.valueless_by_exception());
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

std::vector<double> Weights::to_vector(Type type_) const {
  assert(!_weights.valueless_by_exception());
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

bool Weights::is_constant() const noexcept {
  assert(!_weights.valueless_by_exception());
  return std::holds_alternative<ConstWeight>(_weights);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
bool Weights::is_vector_of_ones() const noexcept {
  assert(!_weights.valueless_by_exception());
  return is_constant() && std::get<ConstWeight>(_weights).w == 1.0;
}

auto Weights::iterator::ConstIt::operator*() const noexcept(ndebug_defined()) -> double {
  assert(!!value);
  bound_check(0, false);
  return value->w;
}

auto Weights::iterator::ConstIt::operator[]([[maybe_unused]] std::ptrdiff_t i_) const
    noexcept(ndebug_defined()) -> double {
  bound_check(i_, false);
  return **this;
}

Weights::iterator::iterator(std::vector<double>::const_iterator it, bool reciprocal)
    : _it(it), _reciprocal(reciprocal) {}

Weights::iterator::iterator(const ConstWeight& weight, std::size_t offset, bool reciprocal)
    : _it(ConstIt{&weight, static_cast<std::ptrdiff_t>(offset)}), _reciprocal(reciprocal) {
  assert(offset <= weight.size);
}

Weights::iterator::iterator(ConstIt it, bool reciprocal) noexcept
    : _it(it), _reciprocal(reciprocal) {}

bool Weights::iterator::is_constant() const noexcept {
  return std::holds_alternative<ConstIt>(_it);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
bool Weights::iterator::operator==(const iterator& other) const noexcept {
  assert(!_it.valueless_by_exception());
  assert(!other._it.valueless_by_exception());
  if (_reciprocal != other._reciprocal) {
    return false;
  }

  if (HICTK_LIKELY(is_constant() == other.is_constant())) {
    return std::visit(
        [&](const auto& it1) {
          using T = remove_cvref_t<decltype(it1)>;
          assert(std::holds_alternative<T>(other._it));
          const auto& it2 = std::get<T>(other._it);
          return it1 == it2;
        },
        _it);
  }
  return false;
}

bool Weights::iterator::operator!=(const iterator& other) const noexcept {
  return !(*this == other);
}

bool Weights::iterator::operator<(const iterator& other) const {
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

bool Weights::iterator::operator<=(const iterator& other) const {
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

bool Weights::iterator::operator>(const iterator& other) const {
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

bool Weights::iterator::operator>=(const iterator& other) const {
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

auto Weights::iterator::operator*() const -> value_type {
  const auto w = std::visit([&](const auto& it) { return *it; }, _it);
  if (_reciprocal) {
    return 1.0 / w;
  }

  return w;
}

auto Weights::iterator::operator[](difference_type i) const -> value_type {
  const auto w = std::visit([&](const auto& it) { return it[i]; }, _it);
  if (_reciprocal) {
    return 1.0 / w;
  }

  return w;
}

auto Weights::iterator::operator++() -> iterator& {
  std::visit([&](auto& it) { ++it; }, _it);
  return *this;
}

auto Weights::iterator::operator++(int) -> iterator {
  return std::visit([&](auto& it_) -> iterator { return {it_++, _reciprocal}; }, _it);
}

auto Weights::iterator::operator+=(difference_type i) -> iterator& {
  std::visit([&](auto& it) { it += i; }, _it);
  return *this;
}

auto Weights::iterator::operator+(difference_type i) const -> iterator {
  return std::visit([&](const auto& it_) -> iterator { return {it_ + i, _reciprocal}; }, _it);
}

auto Weights::iterator::operator--() -> iterator& {
  std::visit([&](auto& it) { --it; }, _it);
  return *this;
}

auto Weights::iterator::operator--(int) -> iterator {
  return std::visit([&](auto& it_) -> iterator { return {it_--, _reciprocal}; }, _it);
}

auto Weights::iterator::operator-=(difference_type i) -> iterator& {
  std::visit([&](auto& it) { it -= i; }, _it);
  return *this;
}

auto Weights::iterator::operator-(difference_type i) const -> iterator {
  return std::visit([&](const auto& it_) -> iterator { return {it_ - i, _reciprocal}; }, _it);
}

auto Weights::iterator::operator-(const iterator& other) const -> difference_type {
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

}  // namespace hictk
