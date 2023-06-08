// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <memory>
#include <type_traits>

#include "hictk/common.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/pixel.hpp"

namespace hictk {

class BinTable;
class Index;

template <typename N, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
class PixelSelector {
  static_assert(std::is_arithmetic_v<N>);

 public:
  class iterator;

 private:
  PixelCoordinates _coord1{};
  PixelCoordinates _coord2{};
  std::shared_ptr<const Index> _index{};
  const Dataset *_pixels_bin1_id{};
  const Dataset *_pixels_bin2_id{};
  const Dataset *_pixels_count{};

 public:
  PixelSelector() = delete;
  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count) noexcept;
  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                PixelCoordinates coords) noexcept;

  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count, PixelCoordinates coord1,
                PixelCoordinates coord2) noexcept;

  template <std::size_t CHUNK_SIZE_OTHER>
  [[nodiscard]] bool operator==(const PixelSelector<N, CHUNK_SIZE_OTHER> &other) const noexcept;
  template <std::size_t CHUNK_SIZE_OTHER>
  [[nodiscard]] bool operator!=(const PixelSelector<N, CHUNK_SIZE_OTHER> &other) const noexcept;

  [[nodiscard]] auto begin() const -> iterator;
  [[nodiscard]] auto end() const -> iterator;

  [[nodiscard]] auto cbegin() const -> iterator;
  [[nodiscard]] auto cend() const -> iterator;

  [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  class iterator {
    using BinIDT = std::uint64_t;
    friend PixelSelector<N, CHUNK_SIZE>;

    Dataset::iterator<BinIDT, CHUNK_SIZE> _bin1_id_it{};
    Dataset::iterator<BinIDT, CHUNK_SIZE> _bin2_id_it{};
    Dataset::iterator<N, CHUNK_SIZE> _count_it{};

    mutable Pixel<N> _value{};
    std::shared_ptr<const Index> _index{};

    PixelCoordinates _coord1{};
    PixelCoordinates _coord2{};

    std::uint64_t _h5_end_offset{};

    explicit iterator(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                      const Dataset &pixels_bin2_id, const Dataset &pixels_count);

    explicit iterator(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                      const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                      PixelCoordinates coord1, PixelCoordinates coord2);

    static auto at_end(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                       const Dataset &pixels_bin2_id, const Dataset &pixels_count) -> iterator;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Pixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;

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

   private:
    void jump_to_row(std::uint64_t bin_id);
    void jump_to_col(std::uint64_t bin_id);
    void jump(std::uint64_t bin1_id, std::uint64_t bin2_id);
    void jump_to_next_overlap();

    [[nodiscard]] std::size_t h5_offset() const noexcept;
    void jump_at_end();
    void refresh();
    void read_pixel() const;

    [[nodiscard]] constexpr bool overlaps_coord1() const noexcept;
    [[nodiscard]] constexpr bool overlaps_coord2() const noexcept;

    [[nodiscard]] bool discard() const;
    constexpr bool is_at_end() const noexcept;
  };
};

}  // namespace hictk

#include "../../../pixel_selector_impl.hpp"
