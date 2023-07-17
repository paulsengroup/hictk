// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5Exception.hpp>
#include <highfive/H5Selection.hpp>
#include <numeric>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/common.hpp"

namespace hictk::cooler {

namespace internal {

// https://www.geeksforgeeks.org/nearest-prime-less-given-number-n/
// https://practice.geeksforgeeks.org/user-profile.php?user=Shashank%20Mishra
template <typename I>
inline bool is_prime(I n) {
  if (n <= 1) {
    return false;
  }
  for (I i = 2; i <= static_cast<I>(std::sqrt(static_cast<double>(n))); ++i) {
    if (n % i == 0) {
      return false;
    }
  }
  return true;
}

template <typename I>
[[nodiscard]] inline I nearest_prime(I n) {
  for (I i = n - 1; i >= 2; --i) {
    if (is_prime(i)) {
      return i;
    }
  }
  return 0;
}

}  // namespace internal

inline HighFive::DataSetCreateProps Dataset::init_create_props(std::uint_fast8_t compression_lvl,
                                                               std::size_t chunk_size) {
  assert(chunk_size != 0);
  HighFive::DataSetCreateProps props{};
  props.add(HighFive::Shuffle());
  props.add(HighFive::Deflate(conditional_static_cast<std::uint32_t>(compression_lvl)));
  props.add(HighFive::Chunking(chunk_size / sizeof(std::int32_t)));
  return props;
}

inline HighFive::DataSetAccessProps Dataset::init_access_props(std::size_t chunk_size,
                                                               std::size_t cache_size, double w0) {
  // https://docs.hdfgroup.org/hdf5/v1_12/group___d_a_p_l.html#ga104d00442c31714ee073dee518f661f
  assert(chunk_size != 0);
  assert(cache_size != 0);

  const auto num_chunks = (std::max)(std::size_t(1), cache_size / chunk_size);
  const auto num_slots = internal::nearest_prime(100 * num_chunks);

  HighFive::DataSetAccessProps props{};
  props.add(HighFive::Caching(num_slots, cache_size, w0));
  return props;
}

inline HighFive::DataSetCreateProps Dataset::default_create_props() {
  return Dataset::init_create_props(DEFAULT_COMPRESSION_LEVEL, DEFAULT_HDF5_CHUNK_SIZE);
}

inline HighFive::DataSetAccessProps Dataset::default_access_props() {
  return Dataset::init_access_props(DEFAULT_HDF5_CHUNK_SIZE, DEFAULT_HDF5_DATASET_CACHE_SIZE,
                                    DEFAULT_HDF5_CACHE_W0);
}

inline Dataset::Dataset(RootGroup root_group, HighFive::DataSet dset)
    : _root_group(std::move(root_group)), _dataset(std::move(dset)) {}

inline Dataset::Dataset(RootGroup root_group, std::string_view path_to_dataset,
                        const HighFive::DataSetAccessProps &aprops)
    : Dataset(root_group, root_group().getDataSet(std::string{path_to_dataset}, aprops)) {}

template <typename T, typename>
inline Dataset::Dataset(RootGroup root_group, std::string_view path_to_dataset,
                        [[maybe_unused]] const T &type, std::size_t max_dim,
                        const HighFive::DataSetAccessProps &aprops,
                        const HighFive::DataSetCreateProps &cprops)
    : Dataset(root_group,
              root_group().createDataSet<T>(std::string{path_to_dataset},
                                            HighFive::DataSpace({0}, {max_dim}), cprops, aprops)) {}

inline Dataset::Dataset(RootGroup root_group, std::string_view path_to_dataset,
                        std::string_view longest_str, std::size_t max_dim,
                        const HighFive::DataSetAccessProps &aprops,
                        const HighFive::DataSetCreateProps &cprops)
    : Dataset(root_group, create_fixed_str_dataset(root_group, path_to_dataset, longest_str.size(),
                                                   max_dim, aprops, cprops)) {}

inline void Dataset::resize(std::size_t new_size) {
  if (new_size > this->_dataset.getElementCount()) {
    this->_dataset.resize({new_size});
  }
}

inline std::pair<std::string, std::string> Dataset::parse_uri(std::string_view uri) {
  const auto pos = uri.rfind('/');
  if (pos == std::string_view::npos) {
    return std::make_pair(std::string{"/"}, std::string{uri});
  }

  if (pos + 1 == uri.size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Invalid dataset URI \"{}\": URI ends with '/'"), uri));
  }

  // clang-format off
  return std::make_pair(std::string{uri.substr(0, pos)},
                        std::string{uri.substr(pos + 1)});
  // clang-format on
}

inline void Dataset::throw_out_of_range_excp(std::size_t offset) const {
  assert(offset >= this->size());

  if (this->empty()) {
    throw std::out_of_range(fmt::format(
        FMT_STRING("Caught an attempt to access an element of dataset {}, which is empty"),
        this->uri(), offset, this->size()));
  }

  throw std::out_of_range(fmt::format(
      FMT_STRING("Caught an attempt to access an element past the end of dataset {} ({} > {})"),
      this->uri(), offset, this->size()));
}

inline void Dataset::throw_out_of_range_excp(std::size_t offset, std::size_t n) const {
  assert(offset + n >= this->size());

  if (this->empty()) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("Caught an attempt to access one or more element(s) of dataset {}, "
                               "which is empty ([{}, {}])"),
                    this->uri(), offset, offset + n));
  }

  throw std::out_of_range(
      fmt::format(FMT_STRING("Caught an attempt to access one or more element(s) past the end of "
                             "dataset {} ([{}-{}] >= {})"),
                  this->uri(), offset, offset + n, this->size()));
}

template <typename T>
inline auto Dataset::make_iterator_at_offset(std::size_t offset, std::size_t chunk_size) const
    -> iterator<T> {
  return iterator<T>(*this, chunk_size, offset);
}

inline HighFive::Selection Dataset::select(std::size_t i) {
  return this->_dataset.select(std::vector<std::size_t>{i});
}

inline HighFive::Selection Dataset::select(std::size_t i) const {
  return this->_dataset.select(std::vector<std::size_t>{i});
}

inline HighFive::Selection Dataset::select(std::size_t i1, std::size_t i2) {
  return this->_dataset.select(std::vector<std::size_t>{i1}, std::vector<std::size_t>{i2});
}

inline HighFive::Selection Dataset::select(std::size_t i1, std::size_t i2) const {
  return this->_dataset.select(std::vector<std::size_t>{i1}, std::vector<std::size_t>{i2});
}

}  // namespace hictk::cooler
