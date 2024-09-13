// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/phmap.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "hictk/pixel.hpp"

namespace hictk::balancing {

class Weights {
 public:
  enum class Type { INFER, DIVISIVE, MULTIPLICATIVE, UNKNOWN };

  class iterator;

 private:
  struct ConstWeight {
    double w{};
    std::size_t size{};
  };
  using WeightVect = std::vector<double>;
  using WeightVectPtr = std::shared_ptr<WeightVect>;

  std::variant<ConstWeight, WeightVectPtr> _weights{ConstWeight{}};
  Type _type{};

 public:
  Weights() = default;
  Weights(std::vector<double> weights, Type type) noexcept;
  Weights(std::vector<double> weights, std::string_view name);
  Weights(double weight, std::size_t size, Type type) noexcept;
  Weights(double weight, std::size_t size, std::string_view name);

  [[nodiscard]] explicit operator bool() const noexcept;

  [[nodiscard]] double at(std::size_t i) const;
  [[nodiscard]] double at(std::size_t i, Type type_) const;

  template <typename N>
  [[nodiscard]] hictk::ThinPixel<N> balance(hictk::ThinPixel<N> p) const;
  template <typename N>
  [[nodiscard]] hictk::Pixel<N> balance(hictk::Pixel<N> p) const;

  [[nodiscard]] Weights operator()(Type type_) const;
  [[nodiscard]] constexpr auto type() const noexcept -> Type;
  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] static auto infer_type(std::string_view name) -> Type;

  void rescale(double scaling_factor) noexcept;
  void rescale(const std::vector<double>& scaling_factors,
               const std::vector<std::uint64_t>& offsets);

  class iterator {
    std::variant<std::vector<double>::const_iterator, const ConstWeight*> _it{nullptr};
    std::ptrdiff_t _i{std::numeric_limits<std::ptrdiff_t>::max()};

    explicit iterator(std::vector<double>::const_iterator it);
    explicit iterator(const ConstWeight& weight);

   public:
  };

 private:
  Weights(std::variant<ConstWeight, WeightVectPtr> weights, Type type_) noexcept;
  template <typename N1, typename N2>
  [[nodiscard]] N1 balance(std::uint64_t bin1_id, std::uint64_t bin2_id, N2 count) const;
  [[nodiscard]] bool is_constant() const noexcept;
};

using WeightMap = phmap::flat_hash_map<std::string, std::shared_ptr<const Weights>>;

}  // namespace hictk::balancing

#include "./impl/weights_impl.hpp"  // NOLINT
