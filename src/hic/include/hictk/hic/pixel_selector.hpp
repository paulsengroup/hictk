// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/hic/block_cache.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/file_reader.hpp"
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
  PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                std::shared_ptr<const internal::HiCFooter> footer_,
                std::shared_ptr<internal::BlockCache> cache_, std::shared_ptr<const BinTable> bins_,
                PixelCoordinates coords) noexcept;

  PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                std::shared_ptr<const internal::HiCFooter> footer_,
                std::shared_ptr<internal::BlockCache> cache_, std::shared_ptr<const BinTable> bins_,
                PixelCoordinates coord1_, PixelCoordinates coord2_) noexcept;

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

  template <typename N>
  [[nodiscard]] std::vector<Pixel<N>> read_all() const;

  [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  [[nodiscard]] MatrixType matrix_type() const noexcept;
  [[nodiscard]] NormalizationMethod normalization() const noexcept;
  [[nodiscard]] MatrixUnit unit() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  [[nodiscard]] const Chromosome &chrom1() const noexcept;
  [[nodiscard]] const Chromosome &chrom2() const noexcept;

  [[nodiscard]] const std::vector<double> &chrom1_norm() const noexcept;
  [[nodiscard]] const std::vector<double> &chrom2_norm() const noexcept;

  [[nodiscard]] const BinTable &bins() const noexcept;
  [[nodiscard]] const internal::HiCFooterMetadata &metadata() const noexcept;

  [[nodiscard]] bool is_inter() const noexcept;
  [[nodiscard]] bool is_intra() const noexcept;
  template <typename N = double>
  [[nodiscard]] N sum() const noexcept;
  [[nodiscard]] double avg() const noexcept;

  [[nodiscard]] std::size_t estimate_optimal_cache_size() const;
  void evict_blocks_from_cache() const;

 private:
  [[nodiscard]] SerializedPixel transform_pixel(SerializedPixel pixel) const;

 public:
  template <typename N>
  class iterator {
    static_assert(std::is_arithmetic_v<N>);
    friend PixelSelector;
    const PixelSelector *_sel{};
    using BufferT = std::vector<Pixel<N>>;
    using BlockIdxBufferT = std::vector<internal::BlockIndex>;

    std::size_t _bin1_id{};
    mutable std::shared_ptr<BlockIdxBufferT> _block_idx_buffer{};
    mutable std::shared_ptr<BufferT> _buffer{};
    mutable std::size_t _buffer_i{};
    mutable std::size_t _pixels_processed{};

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

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;

   private:
    [[nodiscard]] bool is_at_end() const noexcept;
    [[nodiscard]] const BinTable &bins() const noexcept;
    [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
    [[nodiscard]] const PixelCoordinates &coord2() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;

    void read_next_chunk();
    [[nodiscard]] const std::vector<internal::BlockIndex> &find_blocks_overlapping_next_chunk(
        std::size_t num_bins);
    [[nodiscard]] std::size_t compute_chunk_size(double fraction = 0.0005) const noexcept;
  };
};

class PixelSelectorAll {
 public:
  template <typename N>
  class iterator;

 private:
  std::vector<PixelSelector> _selectors{};

 public:
  PixelSelectorAll() = default;
  explicit PixelSelectorAll(std::vector<PixelSelector> selectors_) noexcept;

  template <typename N>
  [[nodiscard]] auto begin() const -> iterator<N>;
  template <typename N>
  [[nodiscard]] auto end() const -> iterator<N>;

  template <typename N>
  [[nodiscard]] auto cbegin() const -> iterator<N>;
  template <typename N>
  [[nodiscard]] auto cend() const -> iterator<N>;

  template <typename N>
  [[nodiscard]] std::vector<Pixel<N>> read_all() const;

  [[nodiscard]] MatrixType matrix_type() const noexcept;
  [[nodiscard]] NormalizationMethod normalization() const noexcept;
  [[nodiscard]] MatrixUnit unit() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] const BinTable &bins() const noexcept;

  template <typename N>
  class iterator {
    const PixelSelectorAll *_sel{};

    using PixelMerger = hictk::internal::PixelMerger<PixelSelector::iterator<N>>;
    std::shared_ptr<PixelMerger> _merger{};
    std::vector<PixelSelector>::const_iterator _it{};
    Pixel<N> _value{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = Pixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    explicit iterator(const PixelSelectorAll &sel);

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;

   private:
    void setup_next_pixel_merger();
  };
};

}  // namespace hictk::hic

#include "../../../pixel_selector_impl.hpp"
