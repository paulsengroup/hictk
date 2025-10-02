// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/pixel_selector.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk {

template <typename PixelSelectorT>
inline PixelSelector::PixelSelector(PixelSelectorT selector,
                                    std::shared_ptr<const balancing::Weights> weights)
    : _sel(std::move(selector)), _weights(std::move(weights)) {
  assert(_weights);
}

template <typename N>
inline auto PixelSelector::begin([[maybe_unused]] bool sorted) const -> iterator<N> {
  return std::visit(
      [&](const auto& sel) {
        using T = std::decay_t<decltype(sel)>;
        if constexpr (std::is_same_v<cooler::PixelSelector, T>) {
          return iterator<N>{sel.template begin<N>(), sel.template end<N>()};
        } else {
          return iterator<N>{sel.template begin<N>(sorted), sel.template end<N>()};
        }
      },
      _sel);
}

template <typename N>
inline auto PixelSelector::end() const -> iterator<N> {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) { return iterator<N>{sel.template end<N>(), sel.template end<N>()}; },
      _sel);
}

template <typename N>
inline auto PixelSelector::cbegin(bool sorted) const -> iterator<N> {
  return begin<N>(sorted);
}

template <typename N>
inline auto PixelSelector::cend() const -> iterator<N> {
  return end<N>();
}

template <typename N>
inline std::vector<Pixel<N>> PixelSelector::read_all() const {
  assert(!_sel.valueless_by_exception());
  return std::visit([&](const auto& sel) { return sel.template read_all<N>(); }, _sel);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline const PixelCoordinates& PixelSelector::coord1() const noexcept {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) -> const PixelCoordinates& {
        using T = std::decay_t<decltype(sel)>;
        if constexpr (std::is_same_v<hic::PixelSelectorAll, T>) {
          static const PixelCoordinates coords{};
          return coords;
        } else {
          return sel.coord1();
        }
      },
      _sel);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline const PixelCoordinates& PixelSelector::coord2() const noexcept {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) -> const PixelCoordinates& {
        using T = std::decay_t<decltype(sel)>;
        if constexpr (std::is_same_v<hic::PixelSelectorAll, T>) {
          static const PixelCoordinates coords{};
          return coords;
        } else {
          return sel.coord2();
        }
      },
      _sel);
}

inline std::uint64_t PixelSelector::size(bool upper_triangle) const {
  assert(!_sel.valueless_by_exception());
  return std::visit([&](const auto& sel) { return sel.size(upper_triangle); }, _sel);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline const BinTable& PixelSelector::bins() const noexcept {
  assert(!_sel.valueless_by_exception());
  return std::visit([&](const auto& sel) -> const BinTable& { return sel.bins(); }, _sel);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
inline std::shared_ptr<const BinTable> PixelSelector::bins_ptr() const noexcept {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) -> std::shared_ptr<const BinTable> { return sel.bins_ptr(); }, _sel);
}

inline PixelSelector PixelSelector::fetch(PixelCoordinates coord1_,
                                          PixelCoordinates coord2_) const {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) -> PixelSelector {
        using T = remove_cvref_t<decltype(sel)>;
        if constexpr (std::is_same_v<hic::PixelSelectorAll, T>) {
          throw std::runtime_error(
              "calling fetch() on a PixelSelector instance set up to fetch genome-wide matrices "
              "from .hic files is not supported");
        } else {
          return PixelSelector{sel.fetch(coord1_, coord2_), _weights};
        }
      },
      _sel);
}

inline const balancing::Weights& PixelSelector::weights() const noexcept {
  assert(_weights);
  return *_weights;
}

template <typename PixelSelectorT>
constexpr const PixelSelectorT& PixelSelector::get() const {
  assert(!_sel.valueless_by_exception());
  return std::get<PixelSelectorT>(_sel);
}

template <typename PixelSelectorT>
constexpr PixelSelectorT& PixelSelector::get() {
  assert(!_sel.valueless_by_exception());
  return std::get<PixelSelectorT>(_sel);
}

constexpr auto PixelSelector::get() const noexcept -> const PixelSelectorVar& { return _sel; }
constexpr auto PixelSelector::get() noexcept -> PixelSelectorVar& { return _sel; }

template <typename N>
template <typename It>
inline PixelSelector::iterator<N>::iterator(It it, It end)
    : _it(std::move(it)), _sentinel(std::move(end)) {}

template <typename N>  // NOLINTNEXTLINE(bugprone-exception-escape)
inline bool PixelSelector::iterator<N>::operator==(const iterator& other) const noexcept {
  assert(!_it.valueless_by_exception());
  return operator_eq(_it, other._it);
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator!=(const iterator& other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator*() const -> const_reference {
  assert(!_it.valueless_by_exception());
  return std::visit([&](const auto& it) -> const_reference { return *it; }, _it);
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator->() const -> const_pointer {
  assert(!_it.valueless_by_exception());
  return std::visit([&](const auto& it) -> const_pointer { return &*it; }, _it);
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++() -> iterator& {
  assert(!_it.valueless_by_exception());
  std::visit([&](auto& it) { ++it; }, _it);
  return *this;
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++(int) -> iterator {
  auto it = *this;
  std::ignore = ++(*this);
  return it;
}

template <typename N>
template <typename IteratorT>
[[nodiscard]] constexpr const IteratorT& PixelSelector::iterator<N>::get() const {
  return std::get<IteratorT>(_it);
}

template <typename N>
template <typename IteratorT>
[[nodiscard]] constexpr IteratorT& PixelSelector::iterator<N>::get() {
  return std::get<IteratorT>(_it);
}

template <typename N>
constexpr auto PixelSelector::iterator<N>::get() const noexcept -> const IteratorVar& {
  return _it;
}

template <typename N>
constexpr auto PixelSelector::iterator<N>::get() noexcept -> IteratorVar& {
  return _it;
}

template <typename N>  // NOLINTNEXTLINE(bugprone-exception-escape)
inline bool PixelSelector::iterator<N>::operator_eq(const IteratorVar& itv1,
                                                    const IteratorVar& itv2) noexcept {
  assert(!itv1.valueless_by_exception());
  assert(!itv2.valueless_by_exception());
  return std::visit(
      [&](const auto& it1) {
        using T = std::decay_t<decltype(it1)>;
        const auto* it2 = std::get_if<T>(&itv2);
        return !!it2 && it1 == *it2;
      },
      itv1);
}

template <typename N>  // NOLINTNEXTLINE(bugprone-exception-escape)
inline bool PixelSelector::iterator<N>::operator_neq(const IteratorVar& itv1,
                                                     const IteratorVar& itv2) noexcept {
  return !(itv1 == itv2);
}

inline File::File(cooler::File clr) : _fp(std::move(clr)) {}
inline File::File(hic::File hf) : _fp(std::move(hf)) {}
inline File::File(std::string_view uri, std::optional<std::uint32_t> resolution_,
                  hic::MatrixType type, hic::MatrixUnit unit) {
  const auto [path, grp] = cooler::parse_cooler_uri(uri);
  if (hic::utils::is_hic_file(path)) {
    *this = File(hic::File(path, resolution_, type, unit));
    return;
  }

  if (type != hic::MatrixType::observed) {
    throw std::runtime_error(
        "matrix type should always be \"observed\" when reading Cooler files.");
  }

  if (unit != hic::MatrixUnit::BP) {
    throw std::runtime_error("matrix unit should always be \"BP\" when reading Cooler files.");
  }

  if (cooler::utils::is_cooler(uri)) {
    *this = File(cooler::File(uri));
    if (resolution_.value_or(resolution()) != resolution()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("found an unexpected resolution while opening file at \"{}\": "
                                 "expected {}, found {}."),
                      uri, *resolution_, resolution()));  // NOLINT(*-unchecked-optional-access)
    }
    return;
  }

  if (!resolution_.has_value()) {
    const auto resolutions = [&]() {
      try {
        return cooler::utils::list_resolutions(uri, false);
      } catch (const std::exception&) {
        return std::vector<std::uint32_t>{};
      }
    }();

    if (resolutions.size() != 1) {
      throw std::runtime_error(
          "resolution is required when opening .mcool files with more than one resolution.");
    }
    resolution_ = resolutions.front();
  }

  const auto new_uri = fmt::format(FMT_STRING("{}::/resolutions/{}"), path, *resolution_);

  if (const auto status = cooler::utils::is_cooler(new_uri); !status) {
    const auto resolution_is_missing =
        status.missing_groups.size() == 1 &&
        status.missing_groups.front().find(
            fmt::format(FMT_STRING("resolutions/{}"), *resolution_)) != std::string::npos;
    if (resolution_is_missing) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("unable to find resolution {} in file \"{}\""), *resolution_, path));
    }
    throw std::runtime_error(
        fmt::format(FMT_STRING("\"{}\" does not look like a valid Cooler file:\n"
                               "Validation report:\n{}"),
                    new_uri, status));
  }

  *this = File(cooler::File(new_uri, cooler::DEFAULT_HDF5_CACHE_SIZE * 4, false));
}

inline std::string File::uri() const {
  assert(!_fp.valueless_by_exception());
  // NOLINTBEGIN(bugprone-branch-clone)
  return std::visit(
      [&](const auto& fp) {
        using T = std::decay_t<decltype(fp)>;
        if constexpr (std::is_same_v<hic::File, T>) {
          return fp.path();
        } else {
          return fp.uri();
        }
      },
      _fp);
  // NOLINTEND(bugprone-branch-clone)
}

inline std::string File::path() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) { return fp.path(); }, _fp);
}

constexpr bool File::is_hic() const noexcept {
  assert(!_fp.valueless_by_exception());
  return std::holds_alternative<hic::File>(_fp);
}

constexpr bool File::is_cooler() const noexcept { return !is_hic(); }

inline auto File::chromosomes() const -> const Reference& {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) -> const Reference& { return fp.chromosomes(); }, _fp);
}

inline auto File::bins() const -> const BinTable& {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) -> const BinTable& { return fp.bins(); }, _fp);
}

inline std::shared_ptr<const BinTable> File::bins_ptr() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& f) -> std::shared_ptr<const BinTable> { return f.bins_ptr(); },
                    _fp);
}

inline std::uint32_t File::resolution() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) { return fp.resolution(); }, _fp);
}

inline std::uint64_t File::nbins() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) { return fp.nbins(); }, _fp);
}

inline std::uint64_t File::nchroms(bool include_ALL) const {
  assert(!_fp.valueless_by_exception());

  return std::visit(
      [&](const auto& fp) {
        using T = remove_cvref_t<decltype(fp)>;
        if constexpr (std::is_same_v<hic::File, T>) {
          return fp.nchroms(include_ALL);
        } else {
          return fp.nchroms();
        }
      },
      _fp);
}

inline PixelSelector File::fetch(const balancing::Method& normalization) const {
  assert(!_fp.valueless_by_exception());
  return std::visit(
      [&](const auto& fp) {
        return PixelSelector{fp.fetch(normalization), fp.normalization_ptr(normalization)};
      },
      _fp);
}

inline PixelSelector File::fetch(std::string_view range, const balancing::Method& normalization,
                                 QUERY_TYPE query_type) const {
  assert(!_fp.valueless_by_exception());
  return fetch(range, range, normalization, query_type);
}

inline PixelSelector File::fetch(std::string_view chrom_name, std::uint32_t start,
                                 std::uint32_t end, const balancing::Method& normalization) const {
  return fetch(chrom_name, start, end, chrom_name, start, end, normalization);
}

inline PixelSelector File::fetch(std::string_view range1, std::string_view range2,
                                 const balancing::Method& normalization,
                                 QUERY_TYPE query_type) const {
  assert(!_fp.valueless_by_exception());
  return std::visit(
      [&](const auto& fp) {
        return PixelSelector{fp.fetch(range1, range2, normalization, query_type),
                             fp.normalization_ptr(normalization)};
      },
      _fp);
}

inline PixelSelector File::fetch(std::string_view chrom1_name, std::uint32_t start1,
                                 std::uint32_t end1, std::string_view chrom2_name,
                                 std::uint32_t start2, std::uint32_t end2,
                                 const balancing::Method& normalization) const {
  assert(!_fp.valueless_by_exception());
  return std::visit(
      [&](const auto& fp) {
        return PixelSelector{
            fp.fetch(chrom1_name, start1, end1, chrom2_name, start2, end2, normalization),
            fp.normalization_ptr(normalization)};
      },
      _fp);
}

inline bool File::has_normalization(std::string_view normalization) const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) { return fp.has_normalization(normalization); }, _fp);
}
inline std::vector<balancing::Method> File::avail_normalizations() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([](const auto& fp) { return fp.avail_normalizations(); }, _fp);
}

inline const balancing::Weights& File::normalization(std::string_view normalization_) const {
  assert(!_fp.valueless_by_exception());
  if (std::holds_alternative<cooler::File>(_fp)) {
    return std::get<cooler::File>(_fp).normalization(normalization_);
  }
  return std::get<hic::File>(_fp).normalization(normalization_);
}

template <typename FileT>
constexpr const FileT& File::get() const {
  return std::get<FileT>(_fp);
}

template <typename FileT>
constexpr FileT& File::get() {
  return std::get<FileT>(_fp);
}

constexpr auto File::get() const noexcept -> const FileVar& { return _fp; }
constexpr auto File::get() noexcept -> FileVar& { return _fp; }

}  // namespace hictk
