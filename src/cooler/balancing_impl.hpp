// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <tsl/hopscotch_map.h>

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace coolerpp {

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

inline Weights::Weights(const BinTable& bins, const Dataset& dset, bool rescale)
    : Weights(bins, dset, Weights::infer_type(dset), rescale) {
  assert(_type != Type::INFER);
  if (_type == Type::UNKNOWN) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), dset.uri()));
  }
}

inline Weights::Weights(const BinTable& bins, const Dataset& dset, Type type, bool rescale)
    : _weights(dset.read_all<std::vector<double>>()), _type(type) {
  if (_type == Type::INFER || type == Type::UNKNOWN) {
    if (dset.has_attribute("divisive_weights")) {
      _type = dset.read_attribute<bool>("divisive_weights") ? Type::DIVISIVE : Type::MULTIPLICATIVE;
    } else {
      _type = Weights::infer_type(dset);
      if (_type == Type::UNKNOWN) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), dset.uri()));
      }
    }
  }

  if (!rescale || !dset.has_attribute("scale")) {
    return;
  }

  const auto cis_only =
      dset.has_attribute("cis_only") ? dset.read_attribute<bool>("cis_only") : false;

  const auto bin_offsets =
      cis_only ? bins.num_bin_prefix_sum() : std::vector<std::uint64_t>{{0, bins.size()}};
  const auto scales = [&]() {
    if (cis_only) {
      std::vector<double> buff;
      dset.read_attribute("scale", buff);
      return buff;
    }
    return std::vector<double>{dset.read_attribute<double>("scale")};
  }();

  assert(!bin_offsets.empty());
  if (bin_offsets.size() - 1 != scales.size()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("failed to read weights from \"{}\": expected {} scale value(s), found {}"),
        dset.uri(), bin_offsets.size() - 1, scales.size()));
  }

  for (std::size_t i = 0; i < scales.size(); ++i) {
    auto first = this->_weights.begin() + std::ptrdiff_t(bin_offsets[i]);
    auto last = this->_weights.begin() + std::ptrdiff_t(bin_offsets[i + 1]);
    std::transform(first, last, first,
                   [s = scales[i]](const double w) { return w * std::sqrt(s); });
  }
}

inline Weights::operator bool() const noexcept { return !this->_weights.empty(); }

inline double Weights::operator[](std::size_t i) const noexcept {
  assert(i < this->_weights.size());
  return this->_weights[i];
}

inline double Weights::at(std::size_t i) const { return this->_weights.at(i); }

inline const std::vector<double>& Weights::operator()() const noexcept { return this->_weights; }

constexpr auto Weights::type() const noexcept -> Type { return this->_type; }

inline auto Weights::infer_type(const Dataset& dset) -> Type {
  auto path = dset.uri();
  auto pos = path.rfind('/');

  if (pos != std::string::npos) {
    path = path.substr(pos + 1);
  }
  return Weights::infer_type(path);
}

inline auto Weights::infer_type(std::string_view name) -> Type {
  const static tsl::hopscotch_map<std::string_view, Type> mappings{
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

template <typename N, std::size_t CHUNK_SIZE>
inline Balancer<N, CHUNK_SIZE>::Balancer(const PixelSelector<N, CHUNK_SIZE>& selector,
                                         std::shared_ptr<const Weights> weights)
    : Balancer(selector.begin(), selector.end(), std::move(weights)) {}

template <typename N, std::size_t CHUNK_SIZE>
inline Balancer<N, CHUNK_SIZE>::Balancer(typename PixelSelector<N, CHUNK_SIZE>::iterator first,
                                         typename PixelSelector<N, CHUNK_SIZE>::iterator last,
                                         std::shared_ptr<const Weights> weights)
    : _first(std::move(first)), _last(std::move(last)), _weights(std::move(weights)) {}

template <typename N, std::size_t CHUNK_SIZE>
inline Weights::Type Balancer<N, CHUNK_SIZE>::type() const noexcept {
  if (!this->_weights) {
    return Weights::Type::UNKNOWN;
  }
  return this->_weights->type();
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto Balancer<N, CHUNK_SIZE>::begin() const -> iterator {
  return this->cbegin();
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto Balancer<N, CHUNK_SIZE>::end() const -> iterator {
  return this->cend();
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto Balancer<N, CHUNK_SIZE>::cbegin() const -> iterator {
  return iterator{this->_first, this->_weights};
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto Balancer<N, CHUNK_SIZE>::cend() const -> iterator {
  return iterator{this->_last, this->_weights};
}

template <typename N, std::size_t CHUNK_SIZE>
inline Balancer<N, CHUNK_SIZE>::iterator::iterator(
    typename PixelSelector<N, CHUNK_SIZE>::iterator it, std::shared_ptr<const Weights> weights)
    : _it(std::move(it)), _weights(std::move(weights)) {}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool Balancer<N, CHUNK_SIZE>::iterator::operator==(
    const Balancer::iterator& other) const noexcept {
  return this->_it == other._it && this->_weights == other._weights;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool Balancer<N, CHUNK_SIZE>::iterator::operator!=(
    const Balancer::iterator& other) const noexcept {
  return !(*this == other);
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool Balancer<N, CHUNK_SIZE>::iterator::operator<(
    const Balancer::iterator& other) const noexcept {
  assert(this->_weights == other._weights);
  return this->_it < other._it;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool Balancer<N, CHUNK_SIZE>::iterator::operator<=(
    const Balancer::iterator& other) const noexcept {
  assert(this->_weights == other._weights);
  return this->_it <= other._it;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool Balancer<N, CHUNK_SIZE>::iterator::operator>(
    const Balancer::iterator& other) const noexcept {
  assert(this->_weights == other._weights);
  return this->_it > other._it;
}

template <typename N, std::size_t CHUNK_SIZE>
constexpr bool Balancer<N, CHUNK_SIZE>::iterator::operator>=(
    const Balancer::iterator& other) const noexcept {
  assert(this->_weights == other._weights);
  return this->_it >= other._it;
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto Balancer<N, CHUNK_SIZE>::iterator::operator*() const -> const_reference {
  const auto& raw_pixel = *this->_it;
  const auto w1 = (*this->_weights)[raw_pixel.coords.bin1.id()];
  const auto w2 = (*this->_weights)[raw_pixel.coords.bin2.id()];
  if (this->_weights->type() == Weights::Type::MULTIPLICATIVE) {
    this->_value = value_type{std::move(raw_pixel.coords),
                              conditional_static_cast<double>(raw_pixel.count) * w1 * w2};
  } else {
    this->_value =
        value_type{std::move(raw_pixel.coords),
                   conditional_static_cast<double>(raw_pixel.count) * (1.0 / w1) * (1.0 / w2)};
  }
  return this->_value;
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto Balancer<N, CHUNK_SIZE>::iterator::operator->() const -> const_pointer {
  return &(*(*this));
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto Balancer<N, CHUNK_SIZE>::iterator::operator++() -> iterator& {
  ++this->_it;

  // signal _value is outdated
  this->_value.count = 0;
  return *this;
}

template <typename N, std::size_t CHUNK_SIZE>
inline auto Balancer<N, CHUNK_SIZE>::iterator::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

}  // namespace coolerpp
