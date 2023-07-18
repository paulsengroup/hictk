// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataType.hpp>

namespace hictk::cooler {

inline HighFive::DataSet Dataset::operator()() { return this->_dataset; }

inline const HighFive::DataSet &Dataset::operator()() const noexcept { return this->_dataset; }

inline std::string Dataset::file_name() const { return this->_root_group().getFile().getName(); }

inline std::string Dataset::hdf5_path() const { return this->_dataset.getPath(); }

inline std::string Dataset::name() const {
  const auto path = hdf5_path();
  const auto last_slash_pos = path.rfind('/');
  if (last_slash_pos == std::string::npos) {
    return path;
  }

  return path.substr(last_slash_pos + 1);
}

inline std::string Dataset::uri() const {
  return fmt::format(FMT_STRING("{}::{}"), this->file_name(), this->hdf5_path());
}

inline std::size_t Dataset::size() const { return this->_dataset.getElementCount(); }

inline bool Dataset::empty() const { return this->size() == 0; }

inline HighFive::DataSet Dataset::get() { return this->_dataset; }
inline const HighFive::DataSet &Dataset::get() const { return this->_dataset; }

inline RootGroup Dataset::get_parent() const { return this->_root_group; }

inline bool Dataset::has_attribute(std::string_view key) const {
  return Attribute::exists(this->_dataset, key);
}

inline HighFive::DataType Dataset::get_h5type() const {
  auto h5type = this->_dataset.getDataType();
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
  auto create_dtype = [&]([[maybe_unused]] auto tunsigned, [[maybe_unused]] auto tsigned) {
    using T1 = decltype(tunsigned);
    using T2 = decltype(tsigned);
    static_assert(std::is_unsigned_v<T1>);
    static_assert(std::is_signed_v<T2>);
    static_assert(sizeof(T1) == sizeof(T2));
    return is_unsigned ? HighFive::create_datatype<T1>() : HighFive::create_datatype<T2>();
  };

  switch (h5type.getSize()) {
    case 1:
      return create_dtype(std::uint8_t{}, std::int8_t{});
    case 2:
      return create_dtype(std::uint16_t{}, std::int16_t{});
    case 4:
      return create_dtype(std::uint32_t{}, std::int32_t{});
    case 8:
      return create_dtype(std::uint64_t{}, std::int64_t{});
    default:
      unreachable_code();
  }
}

template <typename T>
inline auto Dataset::cbegin(std::size_t chunk_size) const -> iterator<T> {
  return iterator<T>(*this, chunk_size);
}

template <typename T>
inline auto Dataset::cend(std::size_t chunk_size) const -> iterator<T> {
  return iterator<T>::make_end_iterator(*this, chunk_size);
}

template <typename T>
inline auto Dataset::begin(std::size_t chunk_size) const -> iterator<T> {
  return this->cbegin<T>(chunk_size);
}

template <typename T>
inline auto Dataset::end(std::size_t chunk_size) const -> iterator<T> {
  return this->cend<T>(chunk_size);
}

}  // namespace hictk::cooler
