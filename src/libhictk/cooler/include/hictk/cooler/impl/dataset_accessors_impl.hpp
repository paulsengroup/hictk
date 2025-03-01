// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#if __has_include(<hdf5/hdf5.h>)
#include <hdf5/H5Ppublic.h>
#include <hdf5/H5Tpublic.h>
#include <hdf5/H5public.h>
#else
#include <H5Ppublic.h>
#include <H5Tpublic.h>
#include <H5public.h>
#endif
#include <fmt/format.h>

#include <cstddef>
#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataType.hpp>
#include <optional>
#include <string>
#include <string_view>

#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"

namespace hictk::cooler {

inline HighFive::DataSet Dataset::operator()() { return _dataset; }

inline const HighFive::DataSet &Dataset::operator()() const noexcept { return _dataset; }

inline std::string Dataset::file_name() const { return _root_group().getFile().getName(); }

inline std::string Dataset::hdf5_path() const { return _dataset.getPath(); }

inline std::string Dataset::name() const {
  auto path = hdf5_path();
  const auto last_slash_pos = path.rfind('/');
  if (last_slash_pos == std::string::npos) {
    return path;
  }

  return path.substr(last_slash_pos + 1);
}

inline std::string Dataset::uri() const {
  return fmt::format(FMT_STRING("{}::{}"), file_name(), hdf5_path());
}

inline std::size_t Dataset::size() const { return _dataset_size; }

inline bool Dataset::empty() const { return size() == 0; }

inline std::size_t Dataset::get_chunk_size() const noexcept { return _chunk_size; }

inline std::size_t Dataset::get_chunk_size(const HighFive::DataSet &dset) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  hsize_t size{};
  const auto dims = H5Pget_chunk(dset.getCreatePropertyList().getId(), 1, &size);
  if (dims != 1) {
    return dset.getElementCount();
  }
  return size;
}

inline HighFive::DataSet Dataset::get() { return _dataset; }
inline const HighFive::DataSet &Dataset::get() const { return _dataset; }

inline RootGroup Dataset::get_parent() const { return _root_group; }

inline bool Dataset::has_attribute(std::string_view key) const {
  return Attribute::exists(_dataset, key);
}

inline HighFive::DataType Dataset::get_h5type() const {
  auto h5type = _dataset.getDataType();
  if (h5type.isFixedLenStr()) {
    return h5type;
  }
  if (h5type.isVariableStr()) {
    return HighFive::create_datatype<std::string>();
  }

  if (h5type.getClass() != HighFive::DataTypeClass::Enum) {
    return h5type;
  }

  // Useful to suppress warnings about treating enum datasets as plain int datasets
  const auto is_unsigned = H5Tget_sign(h5type.getId()) == H5T_SGN_NONE;
  // NOLINTNEXTLINE(*-avoid-non-const-global-variables)
  auto create_dtype = [&]([[maybe_unused]] auto tunsigned, [[maybe_unused]] auto tsigned) {
    using T1 = decltype(tunsigned);
    using T2 = decltype(tsigned);
    static_assert(std::is_unsigned_v<T1>);
    static_assert(std::is_signed_v<T2>);
    static_assert(sizeof(T1) == sizeof(T2));
    return is_unsigned ? HighFive::create_datatype<T1>() : HighFive::create_datatype<T2>();
  };

  switch (h5type.getSize()) {
    case sizeof(std::uint8_t):
      return create_dtype(std::uint8_t{}, std::int8_t{});
    case sizeof(std::uint16_t):
      return create_dtype(std::uint16_t{}, std::int16_t{});
    case sizeof(std::uint32_t):
      return create_dtype(std::uint32_t{}, std::int32_t{});
    case sizeof(std::uint64_t):
      return create_dtype(std::uint64_t{}, std::int64_t{});
    default:
      unreachable_code();
  }
}

template <typename T>
inline auto Dataset::cbegin(std::optional<std::ptrdiff_t> chunk_size) const -> iterator<T> {
  return iterator<T>(*this, chunk_size);
}

template <typename T>
inline auto Dataset::cend(std::optional<std::ptrdiff_t> chunk_size) const -> iterator<T> {
  return iterator<T>::make_end_iterator(*this, chunk_size);
}

template <typename T>
inline auto Dataset::begin(std::optional<std::ptrdiff_t> chunk_size) const -> iterator<T> {
  return cbegin<T>(chunk_size);
}

template <typename T>
inline auto Dataset::end(std::optional<std::ptrdiff_t> chunk_size) const -> iterator<T> {
  return cend<T>(chunk_size);
}

}  // namespace hictk::cooler
