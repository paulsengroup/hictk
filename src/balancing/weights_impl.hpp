// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

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

inline Weights::operator bool() const noexcept { return !this->_weights.empty(); }

inline double Weights::operator[](std::size_t i) const noexcept {
  assert(i < this->_weights.size());
  return this->_weights[i];
}

inline double Weights::at(std::size_t i) const { return this->_weights.at(i); }

template <typename N>
inline ThinPixel<N> Weights::balance(ThinPixel<N> p) const {
  p.count = this->balance<N>(p.bin1_id, p.bin2_id, p.count);
  return p;
}

template <typename N>
inline Pixel<N> Weights::balance(Pixel<N> p) const {
  p.count = this->balance<N>(p.coords.bin1().id(), p.coords.bin2().id(), p.count);
  return p;
}

template <typename N1, typename N2>
inline N1 Weights::balance(std::size_t bin1_id, std::size_t bin2_id, N2 count) const {
  assert(std::is_floating_point_v<N1>);
  const auto w1 = this->_weights[bin1_id];
  const auto w2 = this->_weights[bin2_id];

  auto count_ = conditional_static_cast<double>(count);

  if (this->type() == Weights::Type::MULTIPLICATIVE) {
    count_ *= w1 * w2;
  } else {
    assert(this->type() == Weights::Type::DIVISIVE);
    count_ *= (1.0 / w1) * (1.0 / w2);
  }
  return conditional_static_cast<N1>(count_);
}

inline const std::vector<double>& Weights::operator()() const noexcept { return this->_weights; }

constexpr auto Weights::type() const noexcept -> Type { return this->_type; }

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
       {"weights", Type::MULTIPLICATIVE}}};

  auto it = mappings.find(name);
  if (it == mappings.end()) {
    return Weights::Type::UNKNOWN;
  }
  return it->second;
}

}  // namespace hictk::balancing
