// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <memory>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"

namespace hictk::transformers {

template <typename PixelIt>
class CoarsenPixels {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using N = decltype(std::declval<PixelT>().count);

  PixelIt _first{};
  PixelIt _last{};
  std::shared_ptr<const BinTable> _src_bins{};
  std::shared_ptr<const BinTable> _dest_bins{};
  std::size_t _factor{};

 public:
  class iterator;

  CoarsenPixels(PixelIt first_pixel, PixelIt last_pixel,
                std::shared_ptr<const BinTable> source_bins, std::size_t factor);

  auto begin() const -> iterator;
  auto end() const -> iterator;

  auto cbegin() const -> iterator;
  auto cend() const -> iterator;

  [[nodiscard]] const BinTable &src_bins() const noexcept;
  [[nodiscard]] const BinTable &dest_bins() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> src_bins_ptr() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> dest_bins_ptr() const noexcept;

  [[nodiscard]] auto read_all() const -> std::vector<ThinPixel<N>>;

  class iterator {
    using BufferT = std::vector<ThinPixel<N>>;
    using RowIt = typename BufferT::const_iterator;
    using ColumnMerger = phmap::flat_hash_map<std::uint64_t, BufferT>;

    PixelIt _pixel_it{};
    PixelIt _pixel_last{};
    std::shared_ptr<const BinTable> _src_bins{};
    std::shared_ptr<const BinTable> _dest_bins{};
    std::shared_ptr<BufferT> _buffer{};
    RowIt _it{};

    std::uint64_t _bin1_id_chunk_start{};
    std::uint64_t _bin1_id_chunk_end{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = ThinPixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    iterator(PixelIt first, PixelIt last, std::shared_ptr<const BinTable> src_bins,
             std::shared_ptr<const BinTable> dest_bins);
    static auto at_end(PixelIt last, std::shared_ptr<const BinTable> src_bins,
                       std::shared_ptr<const BinTable> dest_bins) -> iterator;

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    // Current implementation does not allow for an efficient implementation of it++
    auto operator++(int) -> iterator;

   private:
    auto coarsen_chunk_pass1() -> ColumnMerger;
    void coarsen_chunk_pass2(const ColumnMerger &col_merger);
    void process_next_row();
  };
};

}  // namespace hictk::transformers

#include "../../../coarsen_impl.hpp"
