// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// clang-format off
#include "hictk/suppress_warnings.hpp"
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
HICTK_DISABLE_WARNING_POP

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/pixel.hpp"

namespace hictk::balancing {

class Weights {
 public:
  enum class Type : std::uint_fast8_t { INFER, DIVISIVE, MULTIPLICATIVE, UNKNOWN };
  class iterator;

 private:
  struct ConstWeight {
    double w{};
    std::size_t size{};
  };
  using WeightVect = std::vector<double>;
  using WeightVectPtr = std::shared_ptr<WeightVect>;

  std::variant<ConstWeight, WeightVectPtr> _weights{ConstWeight{}};
  Type _type{Type::UNKNOWN};

 public:
  Weights() = default;
  Weights(std::vector<double> weights, Type type);
  Weights(std::vector<double> weights, std::string_view name);
  Weights(double weight, std::size_t size, Type type);
  Weights(double weight, std::size_t size, std::string_view name);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool is_constant() const noexcept;
  [[nodiscard]] bool is_vector_of_ones() const noexcept;

  // NOLINTNEXTLINE(bugprone-exception-escape)
  [[nodiscard]] double operator[](std::size_t i) const noexcept;

  [[nodiscard]] double at(std::size_t i) const;
  [[nodiscard]] double at(std::size_t i, Type type_) const;

  [[nodiscard]] auto begin(Type type = Type::INFER) const -> iterator;
  [[nodiscard]] auto end(Type type = Type::INFER) const -> iterator;

  [[nodiscard]] auto cbegin(Type type = Type::INFER) const -> iterator;
  [[nodiscard]] auto cend(Type type = Type::INFER) const -> iterator;

  template <typename N>
  [[nodiscard]] ThinPixel<N> balance(ThinPixel<N> p) const;
  template <typename N>
  [[nodiscard]] Pixel<N> balance(Pixel<N> p) const;
  template <typename N1, typename N2>
  [[nodiscard]] N1 balance(std::uint64_t bin1_id, std::uint64_t bin2_id, N2 count) const;

  [[nodiscard]] Weights operator()(Type type_) const;
  [[nodiscard]] constexpr auto type() const noexcept -> Type;
  // NOLINTNEXTLINE(bugprone-exception-escape)
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;

  [[nodiscard]] static auto infer_type(std::string_view name) noexcept -> Type;

  void rescale(double scaling_factor);
  void rescale(const std::vector<double> &scaling_factors,
               const std::vector<std::uint64_t> &offsets);

  [[nodiscard]] std::vector<double> to_vector(Type type_ = Type::INFER) const;

  class iterator {
    struct ConstIt {
      const ConstWeight *value{};
      std::ptrdiff_t i{};

      // NOLINTNEXTLINE(bugprone-exception-escape)
      [[nodiscard]] constexpr bool operator==(const ConstIt &other) const noexcept;
      [[nodiscard]] constexpr bool operator!=(const ConstIt &other) const noexcept;
      [[nodiscard]] constexpr bool operator<(const ConstIt &other) const noexcept;
      [[nodiscard]] constexpr bool operator<=(const ConstIt &other) const noexcept;
      [[nodiscard]] constexpr bool operator>(const ConstIt &other) const noexcept;
      [[nodiscard]] constexpr bool operator>=(const ConstIt &other) const noexcept;

      [[nodiscard]] auto operator*() const noexcept(ndebug_defined()) -> double;
      [[nodiscard]] auto operator[](std::ptrdiff_t i_) const noexcept(ndebug_defined()) -> double;

      constexpr auto operator++() noexcept(ndebug_defined()) -> ConstIt &;
      [[nodiscard]] constexpr auto operator++(int) noexcept(ndebug_defined()) -> ConstIt;
      constexpr auto operator+=(std::ptrdiff_t i_) noexcept(ndebug_defined()) -> ConstIt &;
      [[nodiscard]] constexpr auto operator+(std::ptrdiff_t i_) const noexcept(ndebug_defined())
          -> ConstIt;

      constexpr auto operator--() noexcept(ndebug_defined()) -> ConstIt &;
      [[nodiscard]] constexpr auto operator--(int) noexcept(ndebug_defined()) -> ConstIt;
      constexpr auto operator-=(std::ptrdiff_t i_) noexcept(ndebug_defined()) -> ConstIt &;
      [[nodiscard]] constexpr auto operator-(std::ptrdiff_t i_) const noexcept(ndebug_defined())
          -> ConstIt;
      [[nodiscard]] constexpr auto operator-(const ConstIt &other) const noexcept(ndebug_defined())
          -> std::ptrdiff_t;

     private:
      constexpr void bound_check(std::ptrdiff_t offset, bool end_ok = true) const;
    };

    std::variant<std::vector<double>::const_iterator, ConstIt> _it{ConstIt{}};
    bool _reciprocal{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = double;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::random_access_iterator_tag;

    iterator(std::vector<double>::const_iterator it, bool reciprocal);
    iterator(const ConstWeight &weight, std::size_t offset, bool reciprocal);
    iterator(ConstIt it, bool reciprocal) noexcept;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;
    [[nodiscard]] bool operator<(const iterator &other) const;
    [[nodiscard]] bool operator<=(const iterator &other) const;
    [[nodiscard]] bool operator>(const iterator &other) const;
    [[nodiscard]] bool operator>=(const iterator &other) const;

    [[nodiscard]] auto operator*() const -> value_type;
    [[nodiscard]] auto operator[](difference_type i) const -> value_type;

    auto operator++() -> iterator &;
    [[nodiscard]] auto operator++(int) -> iterator;
    auto operator+=(difference_type i) -> iterator &;
    [[nodiscard]] auto operator+(difference_type i) const -> iterator;

    auto operator--() -> iterator &;
    [[nodiscard]] auto operator--(int) -> iterator;
    auto operator-=(difference_type i) -> iterator &;
    auto operator-(difference_type i) const -> iterator;
    [[nodiscard]] auto operator-(const iterator &other) const -> difference_type;

   private:
    [[nodiscard]] bool is_constant() const noexcept;
  };

 private:
  Weights(std::variant<ConstWeight, WeightVectPtr> weights, Type type_) noexcept;
};

using WeightMap = phmap::flat_hash_map<std::string, std::shared_ptr<const Weights>>;

}  // namespace hictk::balancing

#include "./impl/weights_impl.hpp"  // NOLINT
