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

namespace hictk {

inline std::string File::uri() const {
  if (this->hdf5_path() == "/") {
    return this->path();
  }
  return fmt::format(FMT_STRING("{}::{}"), this->path(), this->hdf5_path());
}

inline std::string File::hdf5_path() const { return this->_root_group.hdf5_path(); }
inline std::string File::path() const {
  if (!*this) {
    return "";
  }
  return this->_fp->getName();
}

inline std::uint32_t File::bin_size() const noexcept { return this->_attrs.bin_size; }

inline auto File::chromosomes() const noexcept -> const Reference & {
  return this->bins().chromosomes();
}

inline auto File::bins() const noexcept -> const BinTable & {
  assert(this->_bins);
  return *this->_bins;
}

inline auto File::bins_ptr() const noexcept -> std::shared_ptr<const BinTable> {
  return this->_bins;
}

inline auto File::attributes() const noexcept -> const StandardAttributes & { return this->_attrs; }

inline auto File::group(std::string_view group_name) -> Group & {
  try {
    return this->_groups.at(std::string{group_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Group \"{}\" does not exists!"), group_name));
  }
}

inline auto File::group(std::string_view group_name) const -> const Group & {
  try {
    return this->_groups.at(std::string{group_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Group \"{}\" does not exists!"), group_name));
  }
}

inline auto File::dataset(std::string_view dataset_name) -> Dataset & {
  try {
    if (dataset_name.front() == '/') {
      dataset_name = dataset_name.substr(1);
    }
    return this->_datasets.at(std::string{dataset_name});
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
    return this->_datasets.at(std::string{dataset_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Dataset \"{}\" does not exists!"), dataset_name));
  }
}

inline const internal::NumericVariant &File::pixel_variant() const noexcept {
  return this->_pixel_variant;
}

template <typename T>
inline bool File::has_pixel_of_type() const noexcept {
  return std::holds_alternative<T>(this->_pixel_variant);
}

inline bool File::has_signed_pixels() const noexcept {
  // clang-format off
  return this->has_pixel_of_type<std::int8_t>()  ||
         this->has_pixel_of_type<std::int16_t>() ||
         this->has_pixel_of_type<std::int32_t>() ||
         this->has_pixel_of_type<std::int64_t>();
  // clang-format on
}

inline bool File::has_unsigned_pixels() const noexcept {
  // clang-format off
  return this->has_pixel_of_type<std::uint8_t>()  ||
         this->has_pixel_of_type<std::uint16_t>() ||
         this->has_pixel_of_type<std::uint32_t>() ||
         this->has_pixel_of_type<std::uint64_t>();
  // clang-format on
}

inline bool File::has_integral_pixels() const noexcept {
  return this->has_signed_pixels() || this->has_unsigned_pixels();
}

inline bool File::has_float_pixels() const noexcept {
  // clang-format off
  return this->has_pixel_of_type<float>()  ||
         this->has_pixel_of_type<double>() ||
         this->has_pixel_of_type<long double>();
  // clang-format on
}

template <typename N, std::size_t CHUNK_SIZE>
inline typename PixelSelector<N, CHUNK_SIZE>::iterator File::begin() const {
  // clang-format off
  return PixelSelector<N, CHUNK_SIZE>(this->_index,
                                      this->dataset("pixels/bin1_id"),
                                      this->dataset("pixels/bin2_id"),
                                      this->dataset("pixels/count")
  ).begin();
  // clang-format on
}

template <typename N, std::size_t CHUNK_SIZE>
inline typename PixelSelector<N, CHUNK_SIZE>::iterator File::cbegin() const {
  return this->begin<N, CHUNK_SIZE>();
}

template <typename N, std::size_t CHUNK_SIZE>
inline typename PixelSelector<N, CHUNK_SIZE>::iterator File::end() const {
  // clang-format off
  return PixelSelector<N, CHUNK_SIZE>(this->_index,
                                      this->dataset("pixels/bin1_id"),
                                      this->dataset("pixels/bin2_id"),
                                      this->dataset("pixels/count")
  ).end();
  // clang-format on
}

template <typename N, std::size_t CHUNK_SIZE>
inline typename PixelSelector<N, CHUNK_SIZE>::iterator File::cend() const {
  return this->end<N, CHUNK_SIZE>();
}

inline auto File::index() noexcept -> Index & {
  assert(this->_index);
  return *this->_index;
}

inline auto File::index() const noexcept -> const Index & {
  assert(this->_index);
  return *this->_index;
}

}  // namespace hictk
