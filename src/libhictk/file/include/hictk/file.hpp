// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/pixel_selector.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/tmpdir.hpp"

namespace hictk {

class PixelSelector {
  using PixelSelectorVar =
      std::variant<cooler::PixelSelector, hic::PixelSelector, hic::PixelSelectorAll>;
  PixelSelectorVar _sel{cooler::PixelSelector{}};
  std::shared_ptr<const balancing::Weights> _weights{};

 public:
  template <typename N>
  class iterator;

  PixelSelector() = delete;

  template <typename PixelSelectorT>
  explicit PixelSelector(PixelSelectorT selector,
                         std::shared_ptr<const balancing::Weights> weights);

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

  // NOLINTBEGIN(bugprone-exception-escape)
  [[nodiscard]] const PixelCoordinates &coord1() const noexcept;
  [[nodiscard]] const PixelCoordinates &coord2() const noexcept;

  [[nodiscard]] std::uint64_t size(bool upper_triangle = true) const;

  [[nodiscard]] const BinTable &bins() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;
  // NOLINTEND(bugprone-exception-escape)

  [[nodiscard]] PixelSelector fetch(PixelCoordinates coord1_, PixelCoordinates coord2_) const;

  [[nodiscard]] const balancing::Weights &weights() const noexcept;

  template <typename PixelSelectorT>
  [[nodiscard]] constexpr const PixelSelectorT &get() const;
  template <typename PixelSelectorT>
  [[nodiscard]] constexpr PixelSelectorT &get();
  [[nodiscard]] constexpr auto get() const noexcept -> const PixelSelectorVar &;
  [[nodiscard]] constexpr auto get() noexcept -> PixelSelectorVar &;

  template <typename N>
  class iterator {
    using IteratorVar =
        std::variant<cooler::PixelSelector::iterator<N>, hic::PixelSelector::iterator<N>,
                     hic::PixelSelectorAll::iterator<N>>;
    IteratorVar _it{};
    IteratorVar _sentinel{};

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
    iterator(It it, It end);

    // NOLINTBEGIN(bugprone-exception-escape)
    [[nodiscard]] bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> const_reference;
    [[nodiscard]] auto operator->() const -> const_pointer;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;

    template <typename IteratorT>
    [[nodiscard]] constexpr const IteratorT &get() const;
    template <typename IteratorT>
    [[nodiscard]] constexpr IteratorT &get();
    [[nodiscard]] constexpr auto get() const noexcept -> const IteratorVar &;
    [[nodiscard]] constexpr auto get() noexcept -> IteratorVar &;

   private:
    [[nodiscard]] static bool operator_eq(const IteratorVar &itv1,
                                          const IteratorVar &itv2) noexcept;
    [[nodiscard]] static bool operator_neq(const IteratorVar &itv1,
                                           const IteratorVar &itv2) noexcept;
    // NOLINTEND(bugprone-exception-escape)
  };
};

class File {
  using FileVar = std::variant<cooler::File, hic::File>;
  FileVar _fp{cooler::File{}};

 public:
  using QUERY_TYPE = GenomicInterval::Type;

  explicit File(cooler::File clr);
  explicit File(hic::File hf);
  explicit File(std::string_view uri, std::optional<std::uint32_t> resolution_ = {},
                hic::MatrixType type = hic::MatrixType::observed,
                hic::MatrixUnit unit = hic::MatrixUnit::BP);

  [[nodiscard]] std::string uri() const;
  [[nodiscard]] std::string path() const;

  [[nodiscard]] constexpr bool is_hic() const noexcept;
  [[nodiscard]] constexpr bool is_cooler() const noexcept;

  [[nodiscard]] auto chromosomes() const -> const Reference &;
  [[nodiscard]] auto bins() const -> const BinTable &;
  [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const;

  [[nodiscard]] std::uint32_t resolution() const;
  [[nodiscard]] std::uint64_t nbins() const;
  [[nodiscard]] std::uint64_t nchroms(bool include_ALL = false) const;

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

  [[nodiscard]] bool has_normalization(std::string_view normalization) const;
  [[nodiscard]] std::vector<balancing::Method> avail_normalizations() const;
  [[nodiscard]] const balancing::Weights &normalization(std::string_view normalization_) const;
  [[nodiscard]] std::shared_ptr<const balancing::Weights> normalization_ptr(
      std::string_view normalization_) const;

  template <typename FileT>
  [[nodiscard]] constexpr const FileT &get() const;
  template <typename FileT>
  [[nodiscard]] constexpr FileT &get();
  [[nodiscard]] constexpr auto get() const noexcept -> const FileVar &;
  [[nodiscard]] constexpr auto get() noexcept -> FileVar &;
};

namespace utils {

// NOLINTBEGIN(*-avoid-magic-numbers)

/// Iterable of strings
template <typename N, typename Str>
void merge_to_cool(Str first_uri, Str last_uri, std::string_view dest_uri, std::uint32_t resolution,
                   bool overwrite_if_exists = false, std::size_t chunk_size = 500'000,
                   std::size_t update_frequency = 10'000'000,
                   std::uint32_t compression_lvl = cooler::DEFAULT_COMPRESSION_LEVEL);

/// Iterable of strings
template <typename Str>
void merge_to_hic(
    Str first_file, Str last_file, std::string_view dest_file, std::uint32_t resolution,
    const std::filesystem::path &tmp_dir = internal::TmpDir::default_temp_directory_path(),
    bool overwrite_if_exists = false, std::size_t chunk_size = 500'000, std::size_t n_threads = 1,
    std::uint32_t compression_lvl = 11, bool skip_all_vs_all = false);

// NOLINTEND(*-avoid-magic-numbers)

}  // namespace utils

}  // namespace hictk

#include "./impl/file_impl.hpp"         // NOLINT
#include "./impl/utils_merge_impl.hpp"  // NOLINT
