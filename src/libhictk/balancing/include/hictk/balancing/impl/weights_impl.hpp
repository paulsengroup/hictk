// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::balancing {

inline Weights::Weights(std::vector<double> weights, Type type) noexcept
    : _weights(std::move(weights)), _type(type) {
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

inline Weights::operator bool() const noexcept { return !_weights.empty(); }

inline double Weights::operator[](std::size_t i) const noexcept {
  assert(i < _weights.size());
  return _weights[i];
}

inline double Weights::at(std::size_t i) const { return _weights.at(i); }

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
  const auto w1 = _weights[bin1_id];
  const auto w2 = _weights[bin2_id];

  auto count_ = conditional_static_cast<double>(count);

  if (type() == Weights::Type::MULTIPLICATIVE) {
    count_ *= w1 * w2;
  } else {
    assert(type() == Weights::Type::DIVISIVE);
    count_ *= (1.0 / w1) * (1.0 / w2);
  }
  return conditional_static_cast<N1>(count_);
}

inline const std::vector<double> Weights::operator()(Type type_) const {
  if (type_ != Type::MULTIPLICATIVE && type_ != Type::DIVISIVE) {
    throw std::logic_error("Type should be Type::MULTIPLICATIVE or Type::DIVISIVE");
  }

  if (type_ == type()) {
    return _weights;
  }

  auto weights = _weights;
  std::transform(weights.begin(), weights.end(), weights.begin(),
                 [](const auto n) { return 1.0 / n; });

  return weights;
}

constexpr auto Weights::type() const noexcept -> Type { return _type; }

inline std::size_t Weights::size() const noexcept { return _weights.size(); }

inline auto Weights::infer_type(std::string_view name) -> Type {
  const static phmap::flat_hash_map<std::string_view, Type> mappings{
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

  auto it = mappings.find(name);
  if (it == mappings.end()) {
    return Weights::Type::UNKNOWN;
  }
  return it->second;
}

inline void Weights::rescale(double scaling_factor) noexcept {
  std::transform(_weights.begin(), _weights.end(), _weights.begin(),
                 [&](auto w) { return w * std::sqrt(scaling_factor); });
}

inline void Weights::rescale(const std::vector<double> &scaling_factors,
                             const std::vector<std::uint64_t> &offsets) noexcept {
  for (std::size_t i = 0; i < scaling_factors.size(); ++i) {
    auto first = _weights.begin() + std::ptrdiff_t(offsets[i]);
    auto last = _weights.begin() + std::ptrdiff_t(offsets[i + 1]);
    std::transform(first, last, first,
                   [s = scaling_factors[i]](const double w) { return w * std::sqrt(s); });
  }
}

}  // namespace hictk::balancing
