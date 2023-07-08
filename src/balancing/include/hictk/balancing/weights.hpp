// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/phmap.h>

#include <cstddef>
#include <memory>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::balancing {

class Weights {
 public:
  enum class Type { INFER, DIVISIVE, MULTIPLICATIVE, UNKNOWN };

 private:
  std::vector<double> _weights{};
  Type _type{};

 public:
  Weights() = default;
  Weights(std::vector<double> weights, Type type) noexcept;
  Weights(std::vector<double> weights, std::string_view name);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] double operator[](std::size_t i) const noexcept;

  [[nodiscard]] double at(std::size_t i) const;

  template <typename N>
  [[nodiscard]] hictk::ThinPixel<N> balance(hictk::ThinPixel<N> p) const;
  template <typename N>
  [[nodiscard]] hictk::Pixel<N> balance(hictk::Pixel<N> p) const;

  [[nodiscard]] const std::vector<double>& operator()() const noexcept;
  [[nodiscard]] constexpr auto type() const noexcept -> Type;

  [[nodiscard]] static auto infer_type(std::string_view name) -> Type;

  void rescale(double scaling_factor) noexcept;
  void rescale(const std::vector<double>& slacing_factors,
               const std::vector<std::size_t>& offsets) noexcept;

 private:
  template <typename N1, typename N2>
  [[nodiscard]] N1 balance(std::size_t bin1_id, std::size_t bin2_id, N2 count) const;
};

using WeightMap = phmap::flat_hash_map<std::string, std::shared_ptr<const Weights>>;

}  // namespace hictk::balancing

#include "../../../weights_impl.hpp"
