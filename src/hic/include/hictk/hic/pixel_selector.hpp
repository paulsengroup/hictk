// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "hictk/bin_table.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/hic_file_stream.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic {

class PixelSelector {
  mutable internal::HiCBlockReader _reader{};
  std::shared_ptr<const internal::HiCFooter> _footer{};

  PixelCoordinates _coord1{};
  PixelCoordinates _coord2{};

 public:
  template <typename N>
  class iterator;

  PixelSelector() = delete;
  PixelSelector(std::shared_ptr<internal::HiCFileStream> hfs_,
                std::shared_ptr<const internal::HiCFooter> footer_,
                std::shared_ptr<internal::BlockLRUCache> cache_,
                std::shared_ptr<const BinTable> bins_, PixelCoordinates coords) noexcept;

  PixelSelector(std::shared_ptr<internal::HiCFileStream> hfs_,
                std::shared_ptr<const internal::HiCFooter> footer_,
                std::shared_ptr<internal::BlockLRUCache> cache_,
                std::shared_ptr<const BinTable> bins_, PixelCoordinates coord1_,
                PixelCoordinates coord2_) noexcept;

  [[nodiscard]] bool operator==(const PixelSelector &other) const noexcept;
  [[nodiscard]] bool operator!=(const PixelSelector &other) const noexcept;
  /*
    template <typename N>
    [[nodiscard]] auto begin() const -> iterator<N>;
    template <typename N>
    [[nodiscard]] auto end() const -> iterator<N>;

    template <typename N>
    [[nodiscard]] auto cbegin() const -> iterator<N>;
    template <typename N>
    [[nodiscard]] auto cend() const -> iterator<N>;
    */

  template <typename N>
  [[nodiscard]] std::vector<Pixel<N>> read_all() const;

  [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  [[nodiscard]] const Chromosome &chrom1() const noexcept;
  [[nodiscard]] const Chromosome &chrom2() const noexcept;

  [[nodiscard]] const BinTable &bins() const noexcept;
  [[nodiscard]] const internal::HiCFooterMetadata &metadata() const noexcept;

  [[nodiscard]] bool is_inter() const noexcept;
  [[nodiscard]] bool is_intra() const noexcept;
  template <typename N = double>
  [[nodiscard]] N sum() const noexcept;
  [[nodiscard]] double avg() const noexcept;

 private:
  [[nodiscard]] SerializedPixel process_interaction(SerializedPixel record) const;
  /*
   public:
    template <typename N>
    class iterator {
      static_assert(std::is_arithmetic_v<N>);
      friend PixelSelector;
      const PixelSelector *_sel{};

      std::uint32_t _pos1{};
      std::uint32_t _pos2{};

      std::size_t _col_i{};

      std::shared_ptr<const internal::InteractionBlock> _blk{};
      internal::InteractionBlock::const_iterator _row{};

      mutable Pixel<N> _value{};

     public:
      using difference_type = std::ptrdiff_t;
      using value_type = Pixel<N>;
      using pointer = value_type *;
      using const_pointer = const value_type *;
      using reference = value_type &;
      using const_reference = const value_type &;
      using iterator_category = std::forward_iterator_tag;

      iterator() = default;
      explicit iterator(const PixelSelector &sel);
      [[nodiscard]] static auto at_end(const PixelSelector &sel) -> iterator;

      [[nodiscard]] bool operator==(const iterator &other) const noexcept;
      [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

      [[nodiscard]] bool operator<(const iterator &other) const noexcept;
      [[nodiscard]] bool operator<=(const iterator &other) const noexcept;

      [[nodiscard]] bool operator>(const iterator &other) const noexcept;
      [[nodiscard]] bool operator>=(const iterator &other) const noexcept;

      [[nodiscard]] auto operator*() const -> const_reference;
      // [[nodiscard]] auto operator->() const -> const_pointer;

      auto operator++() -> iterator &;
      auto operator++(int) -> iterator;

     private:
      [[nodiscard]] bool discard() const noexcept;
      [[nodiscard]] bool is_at_end() const noexcept;
      [[nodiscard]] const BinTable &bins() const noexcept;
      [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
      [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

      [[nodiscard]] const internal::InteractionBlock::Row &row() const noexcept;
      [[nodiscard]] std::uint32_t pos1() const noexcept;
      [[nodiscard]] std::uint32_t pos2() const noexcept;
      [[nodiscard]] N count() const noexcept;

      void seek_to_next_block();
      void seek_to_next_overlap() noexcept;
    };
    */
};

}  // namespace hictk::hic

#include "../../../pixel_selector_impl.hpp"
