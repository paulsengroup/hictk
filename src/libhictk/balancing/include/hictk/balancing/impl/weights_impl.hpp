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

inline Weights::Weights(std::vector<double> weights, Type type) noexcept
    : _weights(std::make_shared<WeightVect>(std::move(weights))), _type(type) {
  assert(_type != Type::INFER && _type != Type::UNKNOWN);
}

inline Weights::Weights(std::vector<double> weights, std::string_view name)
    : Weights(std::move(weights), Weights::infer_type(name)) {
  assert(_type != Type::INFER);
  if (_type == Type::UNKNOWN) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), name));
  }
}

inline Weights::Weights(double weight, std::size_t size, Type type) noexcept
    : _weights(ConstWeight{weight, size}), _type(type) {
  assert(_type != Type::INFER && _type != Type::UNKNOWN);
}

inline Weights::Weights(double weight, std::size_t size, std::string_view name)
    : Weights(weight, size, Weights::infer_type(name)) {
  assert(_type != Type::INFER);
  if (_type == Type::UNKNOWN) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), name));
  }
}
inline Weights::Weights(std::variant<ConstWeight, WeightVectPtr> weights, Type type_) noexcept
    : _weights(std::move(weights)), _type(type_) {}

inline Weights::operator bool() const noexcept {
  if (is_constant()) {
    return std::get<ConstWeight>(_weights).size != 0;
  }

  const auto &weights = std::get<WeightVectPtr>(_weights);
  return !!weights && !weights->empty();
}

inline double Weights::at(std::size_t i) const {
  if (is_constant()) {
    const auto &[w, size] = std::get<ConstWeight>(_weights);

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
    count_ *= (1.0 / w1) * (1.0 / w2);
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
    const auto &[w, size] = std::get<ConstWeight>(_weights);
    return {ConstWeight{1.0 / w, size}, type_};
  }

  const auto &buff = std::get<WeightVectPtr>(_weights);
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
  const auto &weights = std::get<WeightVectPtr>(_weights);
  assert(!!weights);
  return weights->size();
}

inline auto Weights::infer_type(std::string_view name) -> Type {
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
                         [&](const auto &p) { return p.first == name; });
  if (it == mappings.end()) {
    return Weights::Type::UNKNOWN;
  }
  return it->second;
}

inline void Weights::rescale(double scaling_factor) noexcept {
  if (is_constant()) {
    auto &w = std::get<ConstWeight>(_weights).w;
    w *= std::sqrt(scaling_factor);
  } else {
    auto &weights = std::get<WeightVectPtr>(_weights);
    assert(!!weights);
    std::transform(weights->begin(), weights->end(), weights->begin(),
                   [&](auto w) { return w * std::sqrt(scaling_factor); });
  }
}

inline void Weights::rescale(const std::vector<double> &scaling_factors,
                             const std::vector<std::uint64_t> &offsets) {
  if (is_constant()) {
    if (scaling_factors.size() == 1) {
      rescale(scaling_factors.front());
      return;
    }
    throw std::runtime_error(
        "rescaling ConstWeight with multiple scaling factors is not supported");
  }

  const auto &weights = std::get<WeightVectPtr>(_weights);
  assert(!!weights);

  for (std::size_t i = 0; i < scaling_factors.size(); ++i) {
    auto first = weights->begin() + std::ptrdiff_t(offsets[i]);
    auto last = weights->begin() + std::ptrdiff_t(offsets[i + 1]);
    std::transform(first, last, first,
                   [s = scaling_factors[i]](const double w) { return w * std::sqrt(s); });
  }
}

inline bool Weights::is_constant() const noexcept {
  return std::holds_alternative<ConstWeight>(_weights);
}

inline Weights::iterator::iterator(std::vector<double>::const_iterator it)
    : _it(std::move(it)), _i(0) {}

inline Weights::iterator::iterator(const hictk::balancing::Weights::ConstWeight &weight)
    : _it(&weight), _i(0) {}

}  // namespace hictk::balancing
