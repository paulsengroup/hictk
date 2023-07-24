// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>

namespace hictk::cooler {

inline std::string File::uri() const {
  if (hdf5_path() == "/") {
    return path();
  }
  return fmt::format(FMT_STRING("{}::{}"), path(), hdf5_path());
}

inline std::string File::hdf5_path() const { return _root_group.hdf5_path(); }
inline std::string File::path() const {
  if (!*this) {
    return "";
  }
  return _root_group().getFile().getName();
}

inline auto File::chromosomes() const noexcept -> const Reference & { return bins().chromosomes(); }

inline auto File::bins() const noexcept -> const BinTable & {
  assert(_bins);
  return *_bins;
}

inline auto File::bins_ptr() const noexcept -> std::shared_ptr<const BinTable> { return _bins; }

inline std::uint32_t File::bin_size() const noexcept { return _attrs.bin_size; }
inline std::uint64_t File::nbins() const { return bins().size(); }
inline std::uint64_t File::nchroms() const { return chromosomes().size(); }
inline std::uint64_t File::nnz() const { return dataset("pixels/count").size(); }

inline auto File::attributes() const noexcept -> const Attributes & { return _attrs; }

inline auto File::group(std::string_view group_name) -> Group & {
  try {
    return _groups.at(std::string{group_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Group \"{}\" does not exists!"), group_name));
  }
}

inline auto File::group(std::string_view group_name) const -> const Group & {
  try {
    return _groups.at(std::string{group_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Group \"{}\" does not exists!"), group_name));
  }
}

inline auto File::dataset(std::string_view dataset_name) -> Dataset & {
  try {
    if (dataset_name.front() == '/') {
      dataset_name = dataset_name.substr(1);
    }
    return _datasets.at(std::string{dataset_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Dataset \"{}\" does not exists!"), dataset_name));
  }
}

inline auto File::dataset(std::string_view dataset_name) const -> const Dataset & {
  try {
    if (dataset_name.front() == '/') {
      dataset_name = dataset_name.substr(1);
    }
    return _datasets.at(std::string{dataset_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Dataset \"{}\" does not exists!"), dataset_name));
  }
}

inline const hictk::internal::NumericVariant &File::pixel_variant() const noexcept {
  return _pixel_variant;
}

template <typename T>
inline bool File::has_pixel_of_type() const noexcept {
  return std::holds_alternative<T>(_pixel_variant);
}

inline bool File::has_signed_pixels() const noexcept {
  // clang-format off
  return has_pixel_of_type<std::int8_t>()  ||
         has_pixel_of_type<std::int16_t>() ||
         has_pixel_of_type<std::int32_t>() ||
         has_pixel_of_type<std::int64_t>();
  // clang-format on
}

inline bool File::has_unsigned_pixels() const noexcept {
  // clang-format off
  return has_pixel_of_type<std::uint8_t>()  ||
         has_pixel_of_type<std::uint16_t>() ||
         has_pixel_of_type<std::uint32_t>() ||
         has_pixel_of_type<std::uint64_t>();
  // clang-format on
}

inline bool File::has_integral_pixels() const noexcept {
  return has_signed_pixels() || has_unsigned_pixels();
}

inline bool File::has_float_pixels() const noexcept {
  // clang-format off
  return has_pixel_of_type<float>()  ||
         has_pixel_of_type<double>() ||
         has_pixel_of_type<long double>();
  // clang-format on
}

template <typename N>
inline PixelSelector::iterator<N> File::begin(std::string_view weight_name) const {
  return fetch(read_weights(balancing::Method(std::string{weight_name}))).template begin<N>();
}

template <typename N>
inline PixelSelector::iterator<N> File::cbegin(std::string_view weight_name) const {
  return begin<N>(weight_name);
}

template <typename N>
inline PixelSelector::iterator<N> File::end(std::string_view weight_name) const {
  return fetch(read_weights(balancing::Method(std::string{weight_name}))).template end<N>();
}

template <typename N>
inline PixelSelector::iterator<N> File::cend(std::string_view weight_name) const {
  return end<N>(weight_name);
}

inline auto File::index() noexcept -> Index & {
  assert(_index);
  return *_index;
}

inline auto File::index() const noexcept -> const Index & {
  assert(_index);
  return *_index;
}

}  // namespace hictk::cooler
