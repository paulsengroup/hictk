// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once
// clang-format off
#include "hictk/suppress_warnings.hpp"
// clang-format on

#include <cstdint>
DISABLE_WARNING_PUSH
DISABLE_WARNING_NULL_DEREF
#include <highfive/H5DataSet.hpp>
DISABLE_WARNING_POP
#include <parallel_hashmap/phmap.h>

#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/generic_variant.hpp"
#include "hictk/variant_buff.hpp"

namespace hictk {

struct RootGroup;

namespace internal {
template <typename T>
struct is_atomic_buffer
    : public std::disjunction<std::is_same<internal::GenericVariant, std::decay_t<T>>,
                              std::is_same<std::string, std::decay_t<T>>,
                              std::is_arithmetic<std::decay_t<T>>> {};

template <typename T>
inline constexpr bool is_atomic_buffer_v = is_atomic_buffer<T>::value;
}  // namespace internal

DISABLE_WARNING_PUSH
DISABLE_WARNING_DEPRECATED_DECLARATIONS
class Dataset {
  RootGroup _root_group{};
  HighFive::DataSet _dataset{};
  mutable internal::VariantBuffer _buff{};

 public:
  template <typename T, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  class iterator;
  template <typename T, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  using const_iterator = iterator<T, CHUNK_SIZE>;

  [[nodiscard]] static HighFive::DataSetCreateProps init_create_props(
      std::uint_fast8_t compression_lvl, std::size_t chunk_size);
  [[nodiscard]] static HighFive::DataSetAccessProps init_access_props(std::size_t chunk_size,
                                                                      std::size_t cache_size,
                                                                      double w0);

  [[nodiscard]] static HighFive::DataSetCreateProps default_create_props();
  [[nodiscard]] static HighFive::DataSetAccessProps default_access_props();

  Dataset() = default;
  Dataset(RootGroup root_group, HighFive::DataSet dset);
  Dataset(RootGroup root_group, std::string_view path_to_dataset,
          const HighFive::DataSetAccessProps &aprops = default_access_props());

  template <typename T, typename = std::enable_if_t<std::is_arithmetic_v<T>>>
  Dataset(RootGroup root_group, std::string_view path_to_dataset, const T &type,
          std::size_t max_dim = HighFive::DataSpace::UNLIMITED,
          const HighFive::DataSetAccessProps &aprops = init_access_props(
              DEFAULT_HDF5_CHUNK_SIZE, DEFAULT_HDF5_CACHE_SIZE * 4, DEFAULT_HDF5_CACHE_W0),
          const HighFive::DataSetCreateProps &cprops = default_create_props());

  Dataset(RootGroup root_group, std::string_view path_to_dataset, std::string_view longest_str,
          std::size_t max_dim = HighFive::DataSpace::UNLIMITED,
          const HighFive::DataSetAccessProps &aprops = init_access_props(
              DEFAULT_HDF5_CHUNK_SIZE, DEFAULT_HDF5_DATASET_CACHE_SIZE * 4, DEFAULT_HDF5_CACHE_W0),
          const HighFive::DataSetCreateProps &cprops = default_create_props());

  const HighFive::DataSet &operator()() const noexcept;
  HighFive::DataSet operator()();

  [[nodiscard]] std::string file_name() const;
  [[nodiscard]] std::string hdf5_path() const;
  [[nodiscard]] std::string uri() const;

  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] bool empty() const;

  [[nodiscard]] HighFive::DataSet get();
  [[nodiscard]] const HighFive::DataSet &get() const;

  [[nodiscard]] RootGroup get_parent() const;

  void resize(std::size_t new_size);

  // Read N values
  template <typename N, typename = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t read(std::vector<N> &buff, std::size_t num, std::size_t offset = 0) const;
  std::size_t read(std::vector<std::string> &buff, std::size_t num, std::size_t offset = 0) const;
  template <std::size_t i = 0>
  std::size_t read(internal::VariantBuffer &vbuff, std::size_t num, std::size_t offset = 0) const;

  template <typename BuffT, typename T = remove_cvref_t<BuffT>,
            typename = std::enable_if_t<!internal::is_atomic_buffer_v<T>>>
  BuffT read_n(std::size_t num, std::size_t offset = 0) const;

  // Read all values
  template <typename BuffT, typename T = remove_cvref_t<BuffT>,
            typename = std::enable_if_t<!internal::is_atomic_buffer_v<T>>>
  std::size_t read_all(BuffT &buff, std::size_t offset = 0) const;

  template <typename BuffT, typename T = remove_cvref_t<BuffT>,
            typename = std::enable_if_t<!internal::is_atomic_buffer_v<T>>>
  BuffT read_all(std::size_t offset = 0) const;

  internal::VariantBuffer read_all(std::size_t offset = 0) const;

  // Read single values
  template <typename N, typename = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t read(N &buff, std::size_t offset) const;
  std::size_t read(std::string &buff, std::size_t offset) const;
  template <std::size_t i = 0>
  std::size_t read(internal::GenericVariant &vbuff, std::size_t offset) const;

  template <typename BuffT, typename T = remove_cvref_t<BuffT>,
            typename = std::enable_if_t<internal::is_atomic_buffer_v<T>>>
  BuffT read(std::size_t offset) const;
  internal::GenericVariant read(std::size_t offset) const;

  template <typename BuffT>
  [[nodiscard]] BuffT read_last() const;
  [[nodiscard]] internal::GenericVariant read_last() const;

  // Write N values
  template <typename N, typename = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t write(const std::vector<N> &buff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);
  std::size_t write(const std::vector<std::string> &buff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);
  std::size_t write(const internal::VariantBuffer &vbuff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);

  template <typename InputIt, typename UnaryOperation = identity,
            typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int>
                * = nullptr>
  std::size_t write(InputIt first_value, InputIt last_value, std::size_t offset = 0,
                    bool allow_dataset_resize = false, UnaryOperation op = identity());

  template <typename InputIt, typename UnaryOperation = identity,
            typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int>
                * = nullptr>
  std::size_t append(InputIt first_value, InputIt last_value, UnaryOperation op = identity());

  // Write single values
  template <typename N, typename = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t write(N buff, std::size_t offset = 0, bool allow_dataset_resize = false);
  std::size_t write(std::string buff, std::size_t offset = 0, bool allow_dataset_resize = false);
  std::size_t write(const internal::GenericVariant &vbuff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);

  template <typename BuffT>
  std::size_t append(const BuffT &buff);

  // Attribute IO
  template <typename T>
  void write_attribute(std::string_view key, const T &value, bool overwrite_if_exists = false);

  template <typename T>
  [[nodiscard]] T read_attribute(std::string_view key) const;

  [[nodiscard]] auto read_attribute(std::string_view key, bool missing_ok = false) const
      -> Attribute::AttributeVar;

  template <typename T>
  void read_attribute(std::string_view key, std::vector<T> &buff) const;

  [[nodiscard]] bool has_attribute(std::string_view key) const;

  template <typename T, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] auto begin() const -> iterator<T, CHUNK_SIZE>;
  template <typename T, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] auto end() const -> iterator<T, CHUNK_SIZE>;

  template <typename T, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] auto cbegin() const -> iterator<T, CHUNK_SIZE>;
  template <typename T, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] auto cend() const -> iterator<T, CHUNK_SIZE>;

  template <typename T, std::size_t CHUNK_SIZE = DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE>
  [[nodiscard]] auto make_iterator_at_offset(std::size_t offset) const -> iterator<T, CHUNK_SIZE>;

  [[nodiscard]] static std::pair<std::string, std::string> parse_uri(std::string_view uri);

 private:
  [[nodiscard]] HighFive::Selection select(std::size_t i);
  [[nodiscard]] HighFive::Selection select(std::size_t i) const;

  [[nodiscard]] HighFive::Selection select(std::size_t i1, std::size_t i2);
  [[nodiscard]] HighFive::Selection select(std::size_t i1, std::size_t i2) const;

  [[nodiscard]] static HighFive::DataSet create_fixed_str_dataset(
      RootGroup &root_grp, std::string_view path, std::size_t max_str_length, std::size_t max_dim,
      const HighFive::DataSetAccessProps &aprops, const HighFive::DataSetCreateProps &cprops);

  [[noreturn]] void throw_out_of_range_excp(std::size_t offset) const;
  [[noreturn]] void throw_out_of_range_excp(std::size_t offset, std::size_t n) const;

  [[nodiscard]] HighFive::DataType get_h5type() const;

 public:
  template <typename T, std::size_t CHUNK_SIZE>
  class iterator {
    static_assert(CHUNK_SIZE != 0);
    friend Dataset;
    mutable std::shared_ptr<std::vector<T>> _buff{};
    const Dataset *_dset{};
    mutable std::size_t _h5_chunk_start{};
    std::size_t _h5_offset{};
#ifndef NDEBUG
    std::size_t _h5_size{};
#endif

    explicit iterator(const Dataset &dset, std::size_t h5_offset = 0, bool init = true);

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::random_access_iterator_tag;

    enum OverlapStatus { UPSTREAM, OVERLAPPING, DOWNSTEAM, UNINITIALIZED };

    iterator() = default;

    [[nodiscard]] constexpr bool operator==(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator!=(const iterator &other) const noexcept;

    [[nodiscard]] constexpr bool operator<(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator<=(const iterator &other) const noexcept;

    [[nodiscard]] constexpr bool operator>(const iterator &other) const noexcept;
    [[nodiscard]] constexpr bool operator>=(const iterator &other) const noexcept;

    [[nodiscard]] auto operator*() const -> value_type;
    [[nodiscard]] auto operator[](std::size_t i) const -> value_type;

    auto operator++() -> iterator &;
    auto operator++(int) -> iterator;
    auto operator+=(std::size_t i) -> iterator &;
    [[nodiscard]] auto operator+(std::size_t i) const -> iterator;

    auto operator--() -> iterator &;
    auto operator--(int) -> iterator;
    auto operator-=(std::size_t i) -> iterator &;
    [[nodiscard]] auto operator-(std::size_t i) const -> iterator;
    [[nodiscard]] auto operator-(const iterator &other) const -> difference_type;

    [[nodiscard]] constexpr std::uint64_t h5_offset() const noexcept;
    [[nodiscard]] constexpr std::size_t underlying_buff_capacity() const noexcept;

    [[nodiscard]] constexpr std::size_t lower_bound() const noexcept;
    [[nodiscard]] constexpr std::size_t upper_bound() const noexcept;

    [[nodiscard]] constexpr auto underlying_buff_status() const noexcept -> OverlapStatus;
    [[nodiscard]] constexpr std::size_t underlying_buff_num_available_rev() const noexcept;
    [[nodiscard]] constexpr std::size_t underlying_buff_num_available_fwd() const noexcept;

    constexpr const Dataset &dataset() const noexcept;

   private:
    void read_chunk_at_offset(std::size_t new_offset) const;

    [[nodiscard]] static constexpr auto make_end_iterator(const Dataset &dset) -> iterator;
  };
};
DISABLE_WARNING_POP

using DatasetMap = phmap::flat_hash_map<std::string, Dataset>;

}  // namespace hictk

#include "../../../dataset_accessors_impl.hpp"
#include "../../../dataset_impl.hpp"
#include "../../../dataset_iterator_impl.hpp"
#include "../../../dataset_read_impl.hpp"
#include "../../../dataset_write_impl.hpp"
