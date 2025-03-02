// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/cooler.hpp"

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <vector>

#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/pixel.hpp"

namespace hictk::cooler {

class PixelSelector {
 public:
  template <typename N>
  class iterator;

 private:
  PixelCoordinates _coord1{};
  PixelCoordinates _coord2{};
  std::shared_ptr<const Index> _index{};
  std::shared_ptr<const BinTable> _bins{};
  const Dataset *_pixels_bin1_id{};
  const Dataset *_pixels_bin2_id{};
  const Dataset *_pixels_count{};
  std::shared_ptr<const balancing::Weights> _weights{};
  bool _symmetric_upper{true};

 public:
  PixelSelector() = default;
  PixelSelector(const Index &index, const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                const Dataset &pixels_count, std::shared_ptr<const balancing::Weights> weights,
                bool symmetric_upper_) noexcept;
  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                const PixelCoordinates &coords, std::shared_ptr<const balancing::Weights> weights,
                bool symmetric_upper_);

  PixelSelector(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                const Dataset &pixels_bin2_id, const Dataset &pixels_count, PixelCoordinates coord1,
                PixelCoordinates coord2, std::shared_ptr<const balancing::Weights> weights,
                bool symmetric_upper_);

  [[nodiscard]] bool operator==(const PixelSelector &other) const noexcept;
  [[nodiscard]] bool operator!=(const PixelSelector &other) const noexcept;

  template <typename N>
  [[nodiscard]] auto begin() const -> iterator<N>;
  template <typename N>
  [[nodiscard]] auto end() const -> iterator<N>;

  template <typename N>
  [[nodiscard]] auto cbegin() const -> iterator<N>;
  template <typename N>
  [[nodiscard]] auto cend() const -> iterator<N>;

  [[nodiscard]] bool empty() const;

  template <typename N>
  [[nodiscard]] std::vector<Pixel<N>> read_all() const;

  [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  [[nodiscard]] const BinTable &bins() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;

  [[nodiscard]] PixelSelector fetch(PixelCoordinates coord1, PixelCoordinates coord2) const;

  [[nodiscard]] const balancing::Weights &weights() const noexcept;

  [[nodiscard]] bool is_symmetric_upper() const noexcept;

  template <typename N>
  class iterator {
    using BinIDT = std::uint64_t;
    friend PixelSelector;

    Dataset::iterator<BinIDT> _bin1_id_it{};
    Dataset::iterator<BinIDT> _bin2_id_it{};
    Dataset::iterator<N> _count_it{};

    mutable ThinPixel<N> _value{};
    std::shared_ptr<const Index> _index{};

    PixelCoordinates _coord1{};
    PixelCoordinates _coord2{};

    std::shared_ptr<const balancing::Weights> _weights{};
    std::uint64_t _h5_end_offset{};
    // this is an offset used to speed up trans queries by skipping over values that are unlikely to
    // overlap with the query range.
    // This optimization works most of the time, and can significantly speed up calls to
    // jump_to_col().
    std::uint64_t _row_head_h5_offset{};

    explicit iterator(const Dataset &pixels_bin1_id, const Dataset &pixels_bin2_id,
                      const Dataset &pixels_count,
                      std::shared_ptr<const balancing::Weights> weights);

    explicit iterator(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                      const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                      PixelCoordinates coord1, PixelCoordinates coord2,
                      std::shared_ptr<const balancing::Weights> weights);

    static auto at_end(std::shared_ptr<const Index> index, const Dataset &pixels_bin1_id,
                       const Dataset &pixels_bin2_id, const Dataset &pixels_count,
                       std::shared_ptr<const balancing::Weights> weights) -> iterator;

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = ThinPixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] bool operator<(const iterator &other) const noexcept;
    [[nodiscard]] bool operator<=(const iterator &other) const noexcept;

    [[nodiscard]] bool operator>(const iterator &other) const noexcept;
    [[nodiscard]] bool operator>=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;

   private:
    void jump_to_row(std::uint64_t bin_id);
    void jump_to_col(std::uint64_t bin_id);
    void jump(std::uint64_t bin1_id, std::uint64_t bin2_id);
    void jump_to_next_row();
    void jump_to_next_overlap();

    [[nodiscard]] std::size_t h5_offset() const noexcept;
    void jump_at_end();
    void refresh();

    [[nodiscard]] constexpr bool overlaps_coord1() const;
    [[nodiscard]] constexpr bool overlaps_coord2() const;

    [[nodiscard]] bool discard() const;
    constexpr bool is_at_end() const;
  };
};

}  // namespace hictk::cooler

#include "./impl/pixel_selector_impl.hpp"  // NOLINT
