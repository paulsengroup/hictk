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
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/pixel_selector.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk {

template <typename PixelSelectorT>
inline PixelSelector::PixelSelector(PixelSelectorT selector) : _sel(std::move(selector)) {}

template <typename N>
inline auto PixelSelector::begin([[maybe_unused]] bool sorted) const -> iterator<N> {
  return std::visit(
      [&](const auto& sel) {
        using T = std::decay_t<decltype(sel)>;
        if constexpr (std::is_same_v<cooler::PixelSelector, T>) {
          return iterator<N>{sel.template begin<N>()};
        } else {
          return iterator<N>{sel.template begin<N>(sorted)};
        }
      },
      _sel);
}

template <typename N>
inline auto PixelSelector::end() const -> iterator<N> {
  return std::visit([&](const auto& sel) { return iterator<N>{sel.template end<N>()}; }, _sel);
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
  return std::visit([&](const auto& sel) { return sel.template read_all<N>(); }, _sel);
}

#ifdef HICTK_WITH_EIGEN
template <typename N>
inline Eigen::SparseMatrix<N> PixelSelector::read_sparse() const {
  return std::visit([&](const auto& sel) { return sel.template read_sparse<N>(); }, _sel);
}

template <typename N>
inline Eigen::Matrix<N, Eigen::Dynamic, Eigen::Dynamic> PixelSelector::read_dense() const {
  return std::visit([&](const auto& sel) { return sel.template read_dense<N>(); }, _sel);
}
#endif

inline const PixelCoordinates& PixelSelector::coord1() const {
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

inline const PixelCoordinates& PixelSelector::coord2() const {
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

inline const BinTable& PixelSelector::bins() const {
  return std::visit([&](const auto& sel) -> const BinTable& { return sel.bins(); }, _sel);
}

template <typename PixelSelectorT>
constexpr const PixelSelectorT& PixelSelector::get() const noexcept {
  return std::get<PixelSelectorT>(_sel);
}

template <typename PixelSelectorT>
constexpr PixelSelectorT& PixelSelector::get() noexcept {
  return std::get<PixelSelectorT>(_sel);
}

constexpr auto PixelSelector::get() const noexcept -> const PixelSelectorVar& { return _sel; }
constexpr auto PixelSelector::get() noexcept -> PixelSelectorVar& { return _sel; }

template <typename N>
template <typename It>
inline PixelSelector::iterator<N>::iterator(It it) : _it(std::move(it)) {}

template <typename N>
inline bool PixelSelector::iterator<N>::operator==(const iterator& other) const noexcept {
  return std::visit(
      [&](const auto& it1) {
        using T = std::decay_t<decltype(it1)>;
        const auto* it2 = std::get_if<T>(&other._it);
        return !!it2 && it1 == *it2;
      },
      _it);
}

template <typename N>
inline bool PixelSelector::iterator<N>::operator!=(const iterator& other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator*() const -> const_reference {
  return std::visit([&](const auto& it) -> const_reference { return *it; }, _it);
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator->() const -> const_pointer {
  return std::visit([&](const auto& it) -> const_pointer { return &*it; }, _it);
}

template <typename N>
inline auto PixelSelector::iterator<N>::operator++() -> iterator& {
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
[[nodiscard]] constexpr const IteratorT& PixelSelector::iterator<N>::get() const noexcept {
  return std::get<IteratorT>(_it);
}

template <typename N>
template <typename IteratorT>
[[nodiscard]] constexpr IteratorT& PixelSelector::iterator<N>::get() noexcept {
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

inline File::File(cooler::File clr) : _fp(std::move(clr)) {}
inline File::File(hic::File hf) : _fp(std::move(hf)) {}
inline File::File(std::string uri, std::uint32_t resolution, hic::MatrixType type,
                  hic::MatrixUnit unit) {
  const auto [path, grp] = cooler::parse_cooler_uri(uri);
  if (hic::utils::is_hic_file(path)) {
    *this = File(hic::File(path, resolution, type, unit));
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
    return;
  }

  *this = File(cooler::File(fmt::format(FMT_STRING("{}::/resolutions/{}"), uri, resolution)));
}

inline std::string File::uri() const {
  return std::visit(
      [&](auto& fp) {
        using T = std::decay_t<decltype(fp)>;
        if constexpr (std::is_same_v<hic::File, T>) {
          return fp.path();
        } else {
          return fp.uri();
        }
      },
      _fp);
}

inline std::string File::path() const {
  return std::visit(
      [&](auto& fp) {
        using T = std::decay_t<decltype(fp)>;
        if constexpr (std::is_same_v<hic::File, T>) {
          return fp.path();
        } else {
          return fp.path();
        }
      },
      _fp);
}

constexpr bool File::is_hic() const noexcept { return std::holds_alternative<hic::File>(_fp); }

constexpr bool File::is_cooler() const noexcept { return !is_hic(); }

inline auto File::chromosomes() const -> const Reference& {
  return std::visit([&](const auto& fp) -> const Reference& { return fp.chromosomes(); }, _fp);
}

inline auto File::bins() const -> const BinTable& {
  return std::visit([&](const auto& fp) -> const BinTable& { return fp.bins(); }, _fp);
}

inline std::shared_ptr<const BinTable> File::bins_ptr() const {
  return std::visit([&](const auto& f) -> std::shared_ptr<const BinTable> { return f.bins_ptr(); },
                    _fp);
}

inline std::uint32_t File::bin_size() const {
  return std::visit([&](const auto& fp) { return fp.bin_size(); }, _fp);
}

inline std::uint64_t File::nbins() const {
  return std::visit([&](const auto& fp) { return fp.nbins(); }, _fp);
}

inline std::uint64_t File::nchroms() const {
  return std::visit([&](const auto& fp) { return fp.nchroms(); }, _fp);
}

inline PixelSelector File::fetch(const balancing::Method& normalization) const {
  return std::visit([&](const auto& fp) { return PixelSelector{fp.fetch(normalization)}; }, _fp);
}

inline PixelSelector File::fetch(std::string_view range, const balancing::Method& normalization,
                                 hictk::File::QUERY_TYPE query_type) const {
  return fetch(range, range, normalization, query_type);
}

inline PixelSelector File::fetch(std::string_view chrom_name, std::uint32_t start,
                                 std::uint32_t end, const balancing::Method& normalization) const {
  return fetch(chrom_name, start, end, chrom_name, start, end, normalization);
}

inline PixelSelector File::fetch(std::string_view range1, std::string_view range2,
                                 const balancing::Method& normalization,
                                 hictk::File::QUERY_TYPE query_type) const {
  return std::visit(
      [&](const auto& fp) {
        return PixelSelector{fp.fetch(range1, range2, normalization, query_type)};
      },
      _fp);
}

inline PixelSelector File::fetch(std::string_view chrom1_name, std::uint32_t start1,
                                 std::uint32_t end1, std::string_view chrom2_name,
                                 std::uint32_t start2, std::uint32_t end2,
                                 const balancing::Method& normalization) const {
  return std::visit(
      [&](const auto& fp) {
        return PixelSelector{
            fp.fetch(chrom1_name, start1, end1, chrom2_name, start2, end2, normalization)};
      },
      _fp);
}

inline bool File::has_normalization(std::string_view normalization) const {
  return std::visit([&](const auto& fp) { return fp.has_normalization(normalization); }, _fp);
}
inline std::vector<balancing::Method> File::avail_normalizations() const {
  return std::visit([](const auto& fp) { return fp.avail_normalizations(); }, _fp);
}

template <typename FileT>
constexpr const FileT& File::get() const noexcept {
  return std::get<FileT>(_fp);
}

template <typename FileT>
constexpr FileT& File::get() noexcept {
  return std::get<FileT>(_fp);
}

constexpr auto File::get() const noexcept -> const FileVar& { return _fp; }
constexpr auto File::get() noexcept -> FileVar& { return _fp; }

}  // namespace hictk
