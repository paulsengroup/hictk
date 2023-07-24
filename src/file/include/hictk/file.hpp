// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <variant>

#include "hictk/cooler/cooler.hpp"
#include "hictk/hic.hpp"

namespace hictk {

class PixelSelector {
  // clang-format off
  std::variant<
      cooler::PixelSelector,
      hic::PixelSelector,
      hic::PixelSelectorAll>
      _sel{cooler::PixelSelector{}};
  // clang-format on

 public:
  template <typename N>
  class iterator;

  template <typename PixelSelectorT>
  explicit PixelSelector(PixelSelectorT selector);

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

  [[nodiscard]] const PixelCoordinates &coord1() const;
  [[nodiscard]] const PixelCoordinates &coord2() const;

  [[nodiscard]] const BinTable &bins() const;

  template <typename PixelSelectorT>
  [[nodiscard]] constexpr PixelSelectorT &get() const noexcept;

  template <typename N>
  class iterator {
    std::variant<cooler::PixelSelector::iterator<N>, hic::PixelSelector::iterator<N>,
                 hic::PixelSelectorAll::iterator<N>>
        _it{};

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = ThinPixel<N>;
    using pointer = value_type *;
    using const_pointer = const value_type *;
    using reference = value_type &;
    using const_reference = const value_type &;
    using iterator_category = std::forward_iterator_tag;

    iterator() = default;
    template <typename It>
    explicit iterator(It it);

    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;

    template <typename IteratorT>
    [[nodiscard]] constexpr IteratorT &get() const noexcept;
  };
};

class File {
  std::variant<cooler::File, hic::File> _fp{cooler::File{}};

 public:
  using QUERY_TYPE = hictk::GenomicInterval::Type;

  explicit File(cooler::File clr);
  explicit File(hic::File hf);
  explicit File(std::string uri, std::uint32_t resolution = 0,
                hic::MatrixType type = hic::MatrixType::observed,
                hic::MatrixUnit unit = hic::MatrixUnit::BP);

  [[nodiscard]] std::string uri() const;
  [[nodiscard]] std::string path() const;

  [[nodiscard]] constexpr bool is_hic() const noexcept ;
  [[nodiscard]] constexpr bool is_cooler() const noexcept ;

  [[nodiscard]] auto chromosomes() const -> const Reference &;
  [[nodiscard]] auto bins() const -> const BinTable &;

  [[nodiscard]] std::uint32_t bin_size() const;
  [[nodiscard]] std::uint64_t nbins() const;
  [[nodiscard]] std::uint64_t nchroms() const;

  [[nodiscard]] PixelSelector fetch(
      const balancing::Method &normalization = balancing::Method::NONE()) const;
  [[nodiscard]] PixelSelector fetch(
      std::string_view range, const balancing::Method &normalization = balancing::Method::NONE(),
      QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  [[nodiscard]] PixelSelector fetch(
      std::string_view chrom_name, std::uint32_t start, std::uint32_t end,
      const balancing::Method &normalization = balancing::Method::NONE()) const;

  [[nodiscard]] PixelSelector fetch(
      std::string_view range1, std::string_view range2,
      const balancing::Method &normalization = balancing::Method::NONE(),
      QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  [[nodiscard]] PixelSelector fetch(
      std::string_view chrom1_name, std::uint32_t start1, std::uint32_t end1,
      std::string_view chrom2_name, std::uint32_t start2, std::uint32_t end2,
      const balancing::Method &normalization = balancing::Method::NONE()) const;

  template <typename FileT>
  [[nodiscard]] constexpr FileT &get() const noexcept;
};

}  // namespace hictk

#include "../../file_impl.hpp"
