// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/dataset.hpp"

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5PropertyList.hpp>
#include <highfive/H5Selection.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace hictk::cooler {

namespace {

// https://www.geeksforgeeks.org/nearest-prime-less-given-number-n/
// https://practice.geeksforgeeks.org/user-profile.php?user=Shashank%20Mishra
template <typename I>
[[nodiscard]] bool is_prime(I n) {
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
[[nodiscard]] I nearest_prime(I n) {
  for (I i = n - 1; i >= 2; --i) {
    if (is_prime(i)) {
      return i;
    }
  }
  return 0;
}

}  // namespace

HighFive::DataSetCreateProps Dataset::init_create_props(std::uint32_t compression_lvl,
                                                        std::size_t chunk_size) {
  assert(chunk_size != 0);
  HighFive::DataSetCreateProps props{};
  props.add(HighFive::Shuffle());
  props.add(HighFive::Deflate(compression_lvl));
  props.add(HighFive::Chunking(chunk_size / sizeof(std::int32_t)));
  return props;
}

HighFive::DataSetAccessProps Dataset::init_access_props(std::size_t chunk_size,
                                                        std::size_t cache_size, double w0) {
  // https://docs.hdfgroup.org/hdf5/v1_12/group___d_a_p_l.html#ga104d00442c31714ee073dee518f661f
  assert(chunk_size != 0);
  assert(cache_size != 0);

  const auto num_chunks = (std::max)(std::size_t{1}, cache_size / chunk_size);
  const auto num_slots = nearest_prime(100 * num_chunks);

  HighFive::DataSetAccessProps props{};
  props.add(HighFive::Caching(num_slots, cache_size, w0));
  return props;
}

HighFive::DataSetCreateProps Dataset::default_create_props() {
  return init_create_props(DEFAULT_COMPRESSION_LEVEL, DEFAULT_HDF5_CHUNK_SIZE);
}

HighFive::DataSetAccessProps Dataset::default_access_props() {
  return init_access_props(DEFAULT_HDF5_CHUNK_SIZE, DEFAULT_HDF5_DATASET_CACHE_SIZE,
                           DEFAULT_HDF5_CACHE_W0);
}

Dataset::Dataset(RootGroup root_group, HighFive::DataSet dset)
    : _root_group(std::move(root_group)),
      _dataset(std::move(dset)),
      _chunk_size(get_chunk_size(_dataset)),
      _dataset_size(_dataset.getElementCount()) {}

Dataset::Dataset(RootGroup root_group, std::string_view path_to_dataset,
                 const HighFive::DataSetAccessProps &aprops)
    : Dataset(root_group, root_group().getDataSet(std::string{path_to_dataset}, aprops)) {}

Dataset::Dataset(RootGroup root_group, std::string_view path_to_dataset,
                 std::string_view longest_str, std::size_t max_dim,
                 const HighFive::DataSetAccessProps &aprops,
                 const HighFive::DataSetCreateProps &cprops)
    : Dataset(root_group, create_fixed_str_dataset(root_group, path_to_dataset, longest_str.size(),
                                                   max_dim, aprops, cprops)) {}

void Dataset::resize(std::size_t new_size) {
  if (new_size > _dataset.getElementCount()) {
    _dataset.resize({new_size});
    _dataset_size = new_size;
  }
}

std::pair<std::string, std::string> Dataset::parse_uri(std::string_view uri) {
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

void Dataset::throw_out_of_range_excp(std::size_t offset) const {
  assert(offset >= size());

  if (empty()) {
    throw std::out_of_range(fmt::format(
        FMT_STRING("Caught an attempt to access an element of dataset {}, which is empty"), uri(),
        offset, size()));
  }

  throw std::out_of_range(fmt::format(
      FMT_STRING("Caught an attempt to access an element past the end of dataset {} ({} > {})"),
      uri(), offset, size()));
}

void Dataset::throw_out_of_range_excp(std::size_t offset, std::size_t n) const {
  assert(offset + n >= size());

  if (empty()) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("Caught an attempt to access one or more element(s) of dataset {}, "
                               "which is empty ([{}, {}])"),
                    uri(), offset, offset + n));
  }

  throw std::out_of_range(
      fmt::format(FMT_STRING("Caught an attempt to access one or more element(s) past the end of "
                             "dataset {} ([{}-{}] >= {})"),
                  uri(), offset, offset + n, size()));
}

const HighFive::Selection &Dataset::select(std::size_t offset, std::size_t count) const {
  _offsets.front() = offset;
  _counts.front() = count;
  _selection.emplace(_dataset.select(_offsets, _counts));
  return *_selection;
}

HighFive::Selection &Dataset::select(std::size_t offset, std::size_t count) {
  _offsets.front() = offset;
  _counts.front() = count;
  _selection.emplace(_dataset.select(_offsets, _counts));
  return *_selection;
}

}  // namespace hictk::cooler
