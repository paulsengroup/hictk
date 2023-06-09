// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/phmap.h>

#include <memory>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/pixel_selector.hpp"

namespace hictk {

class Weights {
 public:
  enum class Type { INFER, DIVISIVE, MULTIPLICATIVE, UNKNOWN };

 private:
  std::vector<double> _weights{};
  Type _type{};

 public:
  Weights() = default;
  explicit Weights(const BinTable &bins, const Dataset &dset, bool rescale = false);
  Weights(const BinTable &bins, const Dataset &dset, Type type, bool rescale = false);
  Weights(std::vector<double> weights, Type type) noexcept;
  Weights(std::vector<double> weights, std::string_view name);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] double operator[](std::size_t i) const noexcept;

  [[nodiscard]] double at(std::size_t i) const;

  [[nodiscard]] const std::vector<double> &operator()() const noexcept;
  [[nodiscard]] constexpr auto type() const noexcept -> Type;

  [[nodiscard]] static auto infer_type(std::string_view name) -> Type;
  [[nodiscard]] static auto infer_type(const Dataset &dset) -> Type;
};

template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
class Balancer {
 public:
  class iterator;

 private:
  typename PixelSelector<N, CHUNK_SIZE>::iterator _first;
  typename PixelSelector<N, CHUNK_SIZE>::iterator _last;
  std::shared_ptr<const Weights> _weights;

 public:
  Balancer() = delete;
  Balancer(const PixelSelector<N, CHUNK_SIZE> &selector, std::shared_ptr<const Weights> weights);
  Balancer(typename PixelSelector<N, CHUNK_SIZE>::iterator first,
           typename PixelSelector<N, CHUNK_SIZE>::iterator last,
           std::shared_ptr<const Weights> weights);

  [[nodiscard]] Weights::Type type() const noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;

  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  class iterator {
    typename PixelSelector<N, CHUNK_SIZE>::iterator _it{};
    std::shared_ptr<const Weights> _weights{};
    mutable Pixel<double> _value{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Pixel<double>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(typename PixelSelector<N, CHUNK_SIZE>::iterator it,
             std::shared_ptr<const Weights> weights);

    [[nodiscard]] constexpr bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] constexpr bool operator<(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator<=(const iterator &other) const noexcept;

    [[nodiscard]] constexpr bool operator>(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator>=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
  };
};

using WeightMap = phmap::flat_hash_map<std::string, std::shared_ptr<const Weights>>;

}  // namespace hictk

#include "../../../balancing_impl.hpp"
