// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5PropertyList.hpp>
#include <highfive/H5Selection.hpp>
#include <iterator>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
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

inline HighFive::DataSetCreateProps Dataset::init_create_props(std::uint32_t compression_lvl,
                                                               std::size_t chunk_size) {
  assert(chunk_size != 0);
  HighFive::DataSetCreateProps props{};
  props.add(HighFive::Shuffle());
  props.add(HighFive::Deflate(compression_lvl));
  props.add(HighFive::Chunking(chunk_size / sizeof(std::int32_t)));
  return props;
}

inline HighFive::DataSetAccessProps Dataset::init_access_props(std::size_t chunk_size,
                                                               std::size_t cache_size, double w0) {
  // https://docs.hdfgroup.org/hdf5/v1_12/group___d_a_p_l.html#ga104d00442c31714ee073dee518f661f
  assert(chunk_size != 0);
  assert(cache_size != 0);

  const auto num_chunks = (std::max)(std::size_t{1}, cache_size / chunk_size);
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
    : _root_group(std::move(root_group)),
      _dataset(std::move(dset)),
      _chunk_size(get_chunk_size(_dataset)),
      _dataset_size(_dataset.getElementCount()) {}

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
  if (new_size > _dataset.getElementCount()) {
    _dataset.resize({new_size});
    _dataset_size = new_size;
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

inline void Dataset::throw_out_of_range_excp(std::size_t offset, std::size_t n) const {
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

template <typename T>
inline auto Dataset::make_iterator_at_offset(std::size_t offset,
                                             std::optional<std::ptrdiff_t> chunk_size) const
    -> iterator<T> {
  return iterator<T>(*this, chunk_size, offset);
}

inline const HighFive::Selection &Dataset::select(std::size_t offset, std::size_t count) const {
  _offsets.front() = offset;
  _counts.front() = count;
  _selection.emplace(_dataset.select(_offsets, _counts));
  return *_selection;
}

inline HighFive::Selection &Dataset::select(std::size_t offset, std::size_t count) {
  _offsets.front() = offset;
  _counts.front() = count;
  _selection.emplace(_dataset.select(_offsets, _counts));
  return *_selection;
}

namespace internal {

template <typename T>
[[nodiscard]] inline std::optional<Dataset::iterator<T>> try_search_in_chunk(
    Dataset::iterator<T> &it, std::size_t i0, std::size_t i1, const T &value) {
  assert(i0 < i1);
  assert(it.buffer().capacity() != 0);

  const auto &chunk = it.buffer();
  assert(i0 >= chunk.start());
  assert(i1 <= chunk.end());

  const auto first_value = chunk[i0];
  const auto last_value = chunk[i1 - 1];

  assert(first_value <= last_value);

  if (first_value <= value && last_value >= value) {
    auto first = it;
    auto last = std::move(it);
    first.seek(i0);
    last.seek(i1);
    return std::lower_bound(std::move(first), std::move(last), value);
  }

  return {};
}

template <typename T>
[[nodiscard]] inline Dataset::iterator<T> compute_pivot(const Dataset::iterator<T> &first,
                                                        const Dataset::iterator<T> &last,
                                                        const T &value,
                                                        bool assume_uniform_distribution) {
  assert(first < last);

  if (!assume_uniform_distribution) {
    // if we can't assume that data is uniformly distributed, simply pick the mid point
    const auto delta = std::distance(first, last);
    return std::next(first, delta / 2);
  }

  const auto first_h5_offset = first.h5_offset();
  const auto last_h5_offset = last.h5_offset();

  const auto first_value = *first;
  const auto last_value = *(last - 1);

  assert(value >= first_value);
  assert(last_value > first_value);

  // compute the ratio between (value - *first) and (*(last - 1) - *first);
  const auto delta = conditional_static_cast<double>(value - first_value);
  const auto range = conditional_static_cast<double>(last_value - first_value);
  const auto cfx = std::clamp(
      conditional_static_cast<double>(delta) / conditional_static_cast<double>(range), 0.0, 1.0);

  // guess a good relative h5_offset to use as pivot.
  // Guessing is done assuming that values between first and last are uniformly distributed
  const auto rel_pivot_offset =
      std::round(cfx * static_cast<double>(last_h5_offset - first_h5_offset));

  const auto pivot_offset = std::clamp(first_h5_offset + static_cast<std::size_t>(rel_pivot_offset),
                                       first_h5_offset + 1, last_h5_offset - 1);
  auto pivot = first;
  pivot.seek(pivot_offset);

  return pivot;
}

template <typename T>  // NOLINTNEXTLINE(readability-function-cognitive-complexity)
inline Dataset::iterator<T> lower_bound_impl(Dataset::iterator<T> first, Dataset::iterator<T> last,
                                             const T &value, bool assume_uniform_distribution,
                                             std::size_t recursion_lvl = 0) {
  assert(first.buffer().capacity() != 0);
  assert(first.buffer().capacity() == last.buffer().capacity());
  assert(first < last);

  const auto first_h5_offset = first.h5_offset();
  const auto last_h5_offset = last.h5_offset();

  const auto &first_chunk = first.buffer();
  const auto &last_chunk = last.buffer();
  assert(!first_chunk.empty());

  if (first_chunk.end() >= last_h5_offset) {
    // only one chunk left

    // we are assigning first to last and then seeking to the appropriate offset just to make sure
    // that first and last have the same underlying buffer.
    // This can sometime improve performance through better cache locality
    last = first;
    last.seek(last_h5_offset);

    return std::lower_bound(std::move(first), std::move(last), value);
  }

  if (*first >= value) {
    // nothing to do: first is downstream of the given value
    return first;
  }

  if (!first_chunk.empty()) {
    // check if value may be contained in the first chunk
    assert(last_h5_offset != 0);
    const auto chunk_h5_last_offset =
        std::clamp(first_chunk.end(), first_chunk.start(), last_h5_offset);
    assert(chunk_h5_last_offset != 0);
    const auto i0 = std::clamp(first_chunk.start(), first_h5_offset, chunk_h5_last_offset);
    const auto i1 = std::clamp(first_chunk.end(), i0, chunk_h5_last_offset);

    auto it = internal::try_search_in_chunk(first, i0, i1, value);
    if (it.has_value()) {
      return *it;
    }
  }

  if (!last_chunk.empty() && last_chunk.start() < last_h5_offset) {
    // check if value may be contained in the last chunk
    assert(last_h5_offset != 0);
    const auto chunk_h5_last_offset =
        std::clamp(last_chunk.end(), last_chunk.start(), last_h5_offset);
    assert(chunk_h5_last_offset != 0);
    const auto i0 = std::clamp(last_chunk.start(), first_h5_offset, chunk_h5_last_offset);
    const auto i1 = std::clamp(last_chunk.end(), i0, chunk_h5_last_offset);

    auto it = internal::try_search_in_chunk(last, i0, i1, value);
    if (it.has_value()) {
      return *it;
    }
  }

  if (recursion_lvl > 4) {
    // if we've not found the correct chunk within the first few attempts, then data is likely not
    // uniformly distributed
    assume_uniform_distribution = false;
  }

  // pick a reasonable pivot point
  auto pivot = internal::compute_pivot(first, last, value, assume_uniform_distribution);

  const auto &pivot_chunk = pivot.buffer();
  assert(!pivot_chunk.empty());
  assert(last_h5_offset != 0);
  assert(pivot_chunk.end() != 0);

  // check if value may be contained in the pivot chunk
  assert(last_h5_offset != 0);
  const auto chunk_h5_last_offset =
      std::clamp(pivot_chunk.end(), pivot_chunk.start(), last_h5_offset - 1);
  assert(chunk_h5_last_offset != 0);
  const auto i0 = std::clamp(pivot_chunk.start(), first_h5_offset, chunk_h5_last_offset);
  const auto i1 = std::clamp(pivot_chunk.end(), i0, chunk_h5_last_offset);

  if (const auto pivot_id = pivot_chunk.id();
      pivot_id != first_chunk.id() && pivot_id != last_chunk.id()) {
    // avoid checking if the given value may be found in the chunk of values that underlie pivot
    // when the chunk is the same as that from the first or last iterator.
    auto it = internal::try_search_in_chunk(pivot, i0, i1, value);
    if (it.has_value()) {
      return *it;
    }
  }

  if (value < pivot_chunk[i0]) {
    // value could be to the left of the pivot point
    pivot.seek(i0);
    return internal::lower_bound_impl(std::move(first), std::move(pivot), value,
                                      assume_uniform_distribution, ++recursion_lvl);
  }

  // value could be to the right of the pivot point
  pivot.seek(i1);
  return internal::lower_bound_impl(std::move(pivot), std::move(last), value,
                                    assume_uniform_distribution, ++recursion_lvl);
}

}  // namespace internal

template <typename T>
inline auto Dataset::lower_bound(iterator<T> first, iterator<T> last, const T &value,
                                 bool assume_uniform_distribution) -> iterator<T> {
  assert(first.dataset().uri() == last.dataset().uri());
  assert(first <= last);

  if (first == last) {
    // nothing to do: range is empty
    return first;
  }

  if (*first >= value) {
    // nothing to do: first is downstream of the given value
    return first;
  }

  // ensure iterators are using the same chunk size
  const auto chunk_size = std::max(first.buffer().capacity(), last.buffer().capacity());
  if (last.buffer().capacity() != chunk_size) {
    last =
        iterator<T>{last._dset, -static_cast<std::ptrdiff_t>(chunk_size), last.h5_offset(), true};
  }

  if (*(last - 1) < value) {
    // value preceding last is upstream of the given value
    return last;
  }

  if (first.buffer().capacity() != chunk_size) {
    first = iterator<T>{first._dset, chunk_size, first.h5_offset(), true};
  }

  return internal::lower_bound_impl(std::move(first), std::move(last), value,
                                    assume_uniform_distribution);
}

}  // namespace hictk::cooler
