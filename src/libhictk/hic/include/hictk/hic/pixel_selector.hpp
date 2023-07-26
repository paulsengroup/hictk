// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICTK_WITH_EIGEN
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#endif
#include <cstddef>
#include <cstdint>
#include <memory>
#include <queue>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/hic/block_reader.hpp"
#include "hictk/hic/cache.hpp"
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

  PixelSelector(const PixelSelector &other) = delete;
  PixelSelector(PixelSelector &&other) = default;

  ~PixelSelector() noexcept;

  PixelSelector &operator=(const PixelSelector &other) = delete;
  PixelSelector &operator=(PixelSelector &&other) = default;

  [[nodiscard]] bool operator==(const PixelSelector &other) const noexcept;
  [[nodiscard]] bool operator!=(const PixelSelector &other) const noexcept;
  template <typename N>
  [[nodiscard]] auto begin(bool sorted = true) const -> iterator<N>;
  template <typename N>
  [[nodiscard]] auto end() const -> iterator<N>;

  template <typename N>
  [[nodiscard]] auto cbegin(bool sorted = true) const -> iterator<N>;
  template <typename N>
  [[nodiscard]] auto cend() const -> iterator<N>;

  template <typename N>
  [[nodiscard]] std::vector<Pixel<N>> read_all() const;

#ifdef HICTK_WITH_EIGEN
  template <typename N>
  [[nodiscard]] Eigen::SparseMatrix<N> read_sparse() const;
  template <typename N>
  [[nodiscard]] Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic> read_dense() const;
#endif

  [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  [[nodiscard]] MatrixType matrix_type() const noexcept;
  [[nodiscard]] balancing::Method normalization() const noexcept;
  [[nodiscard]] MatrixUnit unit() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  [[nodiscard]] const Chromosome &chrom1() const noexcept;
  [[nodiscard]] const Chromosome &chrom2() const noexcept;

  [[nodiscard]] const balancing::Weights &weights1() const noexcept;
  [[nodiscard]] const balancing::Weights &weights2() const noexcept;

  [[nodiscard]] const BinTable &bins() const noexcept;
  [[nodiscard]] const internal::HiCFooterMetadata &metadata() const noexcept;

  [[nodiscard]] bool is_inter() const noexcept;
  [[nodiscard]] bool is_intra() const noexcept;
  template <typename N = double>
  [[nodiscard]] N sum() const noexcept;
  [[nodiscard]] double avg() const noexcept;

  [[nodiscard]] std::size_t estimate_optimal_cache_size(std::size_t num_samples = 500) const;
  void clear_cache() const;

 private:
  template <typename N>
  [[nodiscard]] ThinPixel<N> transform_pixel(ThinPixel<float> pixel) const;

 public:
  template <typename N>
  class iterator {
    static_assert(std::is_arithmetic_v<N>);
    friend PixelSelector;
    using BufferT = std::vector<ThinPixel<N>>;
    using BlockBlacklist = phmap::flat_hash_set<internal::BlockIndex>;

    const PixelSelector *_sel{};
    std::shared_ptr<const internal::Index::Overlap> _block_idx{};
    std::shared_ptr<BlockBlacklist> _block_blacklist{};
    internal::Index::Overlap::const_iterator _block_it{};
    mutable std::shared_ptr<BufferT> _buffer{};
    mutable std::size_t _buffer_i{};
    std::uint32_t _bin1_id{};
    bool _sorted{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = ThinPixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    explicit iterator(const PixelSelector &sel, bool sorted);
    [[nodiscard]] static auto at_end(const PixelSelector &sel) -> iterator;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] bool operator<(const iterator &other) const noexcept;
    [[nodiscard]] bool operator>(const iterator &other) const noexcept;

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
    [[nodiscard]] std::uint64_t bin1_id() const noexcept;
    [[nodiscard]] std::uint64_t bin2_id() const noexcept;

    [[nodiscard]] std::uint32_t compute_chunk_size() const noexcept;
    [[nodiscard]] std::vector<internal::BlockIndex> find_blocks_overlapping_next_chunk(
        std::size_t num_bins);
    void read_next_chunk();
    void read_next_chunk_v9_intra();
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
  [[nodiscard]] auto begin(bool sorted = true) const -> iterator<N>;
  template <typename N>
  [[nodiscard]] auto end() const -> iterator<N>;

  template <typename N>
  [[nodiscard]] auto cbegin(bool sorted = true) const -> iterator<N>;
  template <typename N>
  [[nodiscard]] auto cend() const -> iterator<N>;

  template <typename N>
  [[nodiscard]] std::vector<Pixel<N>> read_all() const;

#ifdef HICTK_WITH_EIGEN
  template <typename N>
  [[nodiscard]] Eigen::SparseMatrix<N> read_sparse() const;
  template <typename N>
  [[nodiscard]] Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic> read_dense() const;
#endif

  [[nodiscard]] MatrixType matrix_type() const noexcept;
  [[nodiscard]] balancing::Method normalization() const noexcept;
  [[nodiscard]] MatrixUnit unit() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] const BinTable &bins() const noexcept;

  template <typename N>
  class iterator {
    struct Pair {
      PixelSelector::iterator<N> first{};  // NOLINT
      PixelSelector::iterator<N> last{};   // NOLINT
      bool operator<(const Pair &other) const noexcept;
      bool operator>(const Pair &other) const noexcept;
    };

    using SelectorQueue = std::queue<const PixelSelector *>;
    using ItPQueue = std::priority_queue<Pair, std::vector<Pair>, std::greater<>>;
    std::shared_ptr<SelectorQueue> _selectors{};
    std::shared_ptr<SelectorQueue> _active_selectors{};
    std::shared_ptr<ItPQueue> _its{};
    bool _sorted{};

    std::uint32_t _chrom1_id{};

    std::shared_ptr<std::vector<ThinPixel<N>>> _buff{};
    std::size_t _i{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = ThinPixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    explicit iterator(const PixelSelectorAll &selector, bool sorted);

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;

   private:
    void init_iterators();
    void read_next_chunk();
  };
};

}  // namespace hictk::hic

#include "./impl/pixel_selector_impl.hpp"
