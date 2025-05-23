// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

// clang-format off
#include "hictk/suppress_warnings.hpp"
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
HICTK_DISABLE_WARNING_POP
// clang-format on

#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <memory>
#include <optional>
#include <queue>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/hic/block_reader.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic {

class PixelSelector {
  std::shared_ptr<internal::HiCBlockReader> _reader{};

  std::shared_ptr<const internal::HiCFooter> _footer{};

  std::shared_ptr<const PixelCoordinates> _coord1{};
  std::shared_ptr<const PixelCoordinates> _coord2{};
  std::optional<std::uint64_t> _diagonal_band_width{};

 public:
  template <typename N>
  class iterator;

  PixelSelector() = default;
  PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                std::shared_ptr<const internal::HiCFooter> footer_,
                std::shared_ptr<internal::BlockCache> cache_, std::shared_ptr<const BinTable> bins_,
                const PixelCoordinates &coords, std::optional<std::uint64_t> diagonal_band_width);

  PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                std::shared_ptr<const internal::HiCFooter> footer_,
                std::shared_ptr<internal::BlockCache> cache_, std::shared_ptr<const BinTable> bins_,
                PixelCoordinates coord1_, PixelCoordinates coord2_,
                std::optional<std::uint64_t> diagonal_band_width);

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

  [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  [[nodiscard]] std::uint64_t size(bool upper_triangle = true) const;

  [[nodiscard]] MatrixType matrix_type() const noexcept;
  [[nodiscard]] const balancing::Method &normalization() const noexcept;
  [[nodiscard]] MatrixUnit unit() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  [[nodiscard]] const Chromosome &chrom1() const noexcept;
  [[nodiscard]] const Chromosome &chrom2() const noexcept;

  [[nodiscard]] const balancing::Weights &weights1() const noexcept;
  [[nodiscard]] const balancing::Weights &weights2() const noexcept;

  [[nodiscard]] const BinTable &bins() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;
  [[nodiscard]] const internal::HiCFooterMetadata &metadata() const noexcept;

  [[nodiscard]] bool is_inter() const noexcept;
  [[nodiscard]] bool is_intra() const noexcept;
  [[nodiscard]] bool empty() const noexcept;

  // NOLINTNEXTLINE(*-avoid-magic-numbers)
  [[nodiscard]] std::size_t estimate_optimal_cache_size(std::size_t num_samples = 500) const;
  void clear_cache() const;

  [[nodiscard]] PixelSelector fetch(PixelCoordinates coord1_, PixelCoordinates coord2_) const;

 private:
  PixelSelector(std::shared_ptr<internal::HiCBlockReader> reader_,
                std::shared_ptr<const internal::HiCFooter> footer_,
                std::shared_ptr<const PixelCoordinates> coord1_,
                std::shared_ptr<const PixelCoordinates> coord2_,
                std::optional<std::uint64_t> diagonal_band_width);
  template <typename N>
  [[nodiscard]] ThinPixel<N> transform_pixel(ThinPixel<float> pixel) const;

 public:
  template <typename N>
  class iterator {
    static_assert(std::is_arithmetic_v<N>);
    friend PixelSelector;
    using BufferT = std::vector<ThinPixel<N>>;
    using BlockBlacklist = phmap::flat_hash_set<internal::BlockIndex>;

    std::shared_ptr<internal::HiCBlockReader> _reader{};
    std::shared_ptr<const PixelCoordinates> _coord1{};
    std::shared_ptr<const PixelCoordinates> _coord2{};
    std::shared_ptr<const internal::HiCFooter> _footer{};
    std::shared_ptr<const internal::Index::Overlap> _block_idx{};
    std::shared_ptr<BlockBlacklist> _block_blacklist{};
    internal::Index::Overlap::const_iterator _block_it{};
    mutable std::shared_ptr<BufferT> _buffer{};
    mutable std::size_t _buffer_i{};
    std::uint32_t _bin1_id{};
    std::uint64_t _diagonal_band_width{};
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
    explicit iterator(const PixelSelector &sel, bool sorted,
                      std::optional<std::uint64_t> diagonal_band_width);
    [[nodiscard]] static auto at_end(std::shared_ptr<internal::HiCBlockReader> reader,
                                     std::shared_ptr<const PixelCoordinates> coord1,
                                     std::shared_ptr<const PixelCoordinates> coord2) -> iterator;

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
        std::size_t num_bins) const;

    void read_next_chunk();
    void read_next_chunk_unsorted();
    void read_next_chunk_sorted();
    void read_next_chunk_v9_intra_sorted();
    [[nodiscard]] ThinPixel<N> transform_pixel(ThinPixel<float> pixel) const;
    [[nodiscard]] static std::shared_ptr<const internal::Index::Overlap> preload_block_index(
        const PixelSelector &sel, std::optional<std::uint64_t> diagonal_band_width, bool sorted);
  };
};

class PixelSelectorAll {
 public:
  template <typename N>
  class iterator;

 private:
  std::vector<PixelSelector> _selectors{};
  std::shared_ptr<const BinTable> _bins{};
  mutable std::shared_ptr<internal::WeightCache> _weight_cache{};

 public:
  PixelSelectorAll() = default;
  explicit PixelSelectorAll(std::vector<PixelSelector> selectors_,
                            std::shared_ptr<internal::WeightCache> weight_cache = nullptr);
  explicit PixelSelectorAll(std::shared_ptr<const BinTable> bins_,
                            std::shared_ptr<internal::WeightCache> weight_cache = nullptr) noexcept;

  [[nodiscard]] bool empty() const noexcept;

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

  [[nodiscard]] std::uint64_t size(bool upper_triangle = true) const noexcept;

  [[nodiscard]] MatrixType matrix_type() const noexcept;
  [[nodiscard]] const balancing::Method &normalization() const noexcept;
  [[nodiscard]] MatrixUnit unit() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;
  [[nodiscard]] const BinTable &bins() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;
  [[nodiscard]] const balancing::Weights &weights() const;

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
    iterator(const iterator &other);
    iterator(iterator &&other) noexcept = default;

    ~iterator() noexcept = default;

    auto operator=(const iterator &other) -> iterator &;
    auto operator=(iterator &&other) noexcept -> iterator & = default;

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

#include "./impl/pixel_selector_impl.hpp"  // NOLINT
