// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/cooler.hpp"

// clang-format off
#include "hictk/suppress_warnings.hpp"
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
HICTK_DISABLE_WARNING_POP

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_NULL_DEREFERENCE
#include <highfive/H5DataSet.hpp>
HICTK_DISABLE_WARNING_POP
// clang-format on

#include <cstddef>
#include <cstdint>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5DataType.hpp>
#include <highfive/H5PropertyList.hpp>
#include <highfive/H5Selection.hpp>
#include <iterator>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/common.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/generic_variant.hpp"
#include "hictk/type_traits.hpp"
#include "hictk/variant_buff.hpp"

namespace hictk::cooler {
namespace internal {
template <typename T>
struct is_scalar_buffer
    : std::disjunction<std::is_same<hictk::internal::GenericVariant, std::decay_t<T>>,
                       std::is_same<std::string, std::decay_t<T>>,
                       std::is_arithmetic<std::decay_t<T>>> {};

template <typename T>
inline constexpr bool is_scalar_buffer_v = is_scalar_buffer<T>::value;

template <typename T>
class COWChunk {
  using BufferT = std::vector<T>;
  using SharedBufferT = std::shared_ptr<BufferT>;
  SharedBufferT _buff{};
  std::size_t _start{};

  static inline const std::vector<T> _empty_buffer{};

 public:
  COWChunk() noexcept = default;
  COWChunk(std::size_t start_, SharedBufferT data_, std::size_t capacity_ = 0) noexcept;
  COWChunk(std::size_t start_, BufferT data_, std::size_t capacity_ = 0);

  [[nodiscard]] constexpr std::size_t id() const noexcept;
  [[nodiscard]] constexpr std::size_t start() const noexcept;
  [[nodiscard]] std::size_t end() const noexcept;

  [[nodiscard]] std::size_t capacity() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t use_count() const noexcept;

  [[nodiscard]] auto operator()() const noexcept -> const BufferT &;
  [[nodiscard]] auto operator()() noexcept -> BufferT &;
  // The indices refer to the whole sequence, not just values in the chunk itself
  [[nodiscard]] auto operator()(std::size_t i) const noexcept -> std::optional<T>;
  [[nodiscard]] auto operator[](std::size_t i) const noexcept -> T;

  void update(std::size_t start_) noexcept;
  void update(std::size_t start_, SharedBufferT data_);
  void update(std::size_t start_, BufferT data_);
  void resize(std::size_t new_size, bool shrink_to_fit = false);
  void reserve(std::size_t new_capacity);
  void reset_buffer() noexcept;
};

}  // namespace internal

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
class Dataset {
  using VariantBuffer = hictk::internal::VariantBuffer;
  using GenericVariant = hictk::internal::GenericVariant;
  RootGroup _root_group{};
  HighFive::DataSet _dataset{};
  mutable std::vector<std::size_t> _offsets{0};
  mutable std::vector<std::size_t> _counts{0};
  mutable std::optional<HighFive::Selection> _selection{};
  mutable VariantBuffer _buff{};
  std::size_t _chunk_size{};
  std::size_t _dataset_size{};

 public:
  template <typename T>
  class iterator;
  template <typename T>
  using const_iterator = iterator<T>;

  [[nodiscard]] static HighFive::DataSetCreateProps init_create_props(std::uint32_t compression_lvl,
                                                                      std::size_t chunk_size_);
  [[nodiscard]] static HighFive::DataSetAccessProps init_access_props(std::size_t chunk_size_,
                                                                      std::size_t cache_size,
                                                                      double w0);

  [[nodiscard]] static HighFive::DataSetCreateProps default_create_props();
  [[nodiscard]] static HighFive::DataSetAccessProps default_access_props();

  Dataset() = default;
  Dataset(RootGroup root_group, HighFive::DataSet dset);
  Dataset(RootGroup root_group, std::string_view path_to_dataset,
          const HighFive::DataSetAccessProps &aprops = default_access_props());

  // NOLINTNEXTLINE(*-avoid-c-arrays)
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
  [[nodiscard]] std::string name() const;
  [[nodiscard]] std::string uri() const;

  [[nodiscard]] std::size_t size() const;
  [[nodiscard]] bool empty() const;

  [[nodiscard]] std::size_t get_chunk_size() const noexcept;

  [[nodiscard]] HighFive::DataSet get();
  [[nodiscard]] const HighFive::DataSet &get() const;

  [[nodiscard]] RootGroup get_parent() const;

  void resize(std::size_t new_size);

  // Read N values
  template <typename N, typename = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t read(std::vector<N> &buff, std::size_t num, std::size_t offset = 0) const;
  std::size_t read(std::vector<std::string> &buff, std::size_t num, std::size_t offset = 0) const;
  template <std::size_t i = 0>
  std::size_t read(VariantBuffer &vbuff, std::size_t num, std::size_t offset = 0) const;

  template <typename BuffT, typename T = remove_cvref_t<BuffT>,
            typename = std::enable_if_t<!internal::is_scalar_buffer_v<T>>>
  BuffT read_n(std::size_t num, std::size_t offset = 0) const;

  // Read all values
  template <typename BuffT, typename T = remove_cvref_t<BuffT>,
            typename = std::enable_if_t<!internal::is_scalar_buffer_v<T>>>
  std::size_t read_all(BuffT &buff, std::size_t offset = 0) const;

  template <typename BuffT, typename T = remove_cvref_t<BuffT>,
            typename = std::enable_if_t<!internal::is_scalar_buffer_v<T>>>
  BuffT read_all(std::size_t offset = 0) const;

  VariantBuffer read_all(std::size_t offset = 0) const;

  // Read single values
  template <typename N, typename = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t read(N &buff, std::size_t offset) const;
  std::size_t read(std::string &buff, std::size_t offset) const;
  template <std::size_t i = 0>
  std::size_t read(GenericVariant &vbuff, std::size_t offset) const;

  template <typename BuffT, typename T = remove_cvref_t<BuffT>,
            typename = std::enable_if_t<internal::is_scalar_buffer_v<T>>>
  BuffT read(std::size_t offset) const;
  GenericVariant read(std::size_t offset) const;

  template <typename BuffT>
  [[nodiscard]] BuffT read_last() const;
  [[nodiscard]] GenericVariant read_last() const;

  // Write N values
  template <typename N, typename = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t write(const std::vector<N> &buff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);
  std::size_t write(const std::vector<std::string> &buff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);
  std::size_t write(const VariantBuffer &vbuff, std::size_t offset = 0,
                    bool allow_dataset_resize = false);

  template <typename InputIt, typename UnaryOperation = identity,
            typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int>
                * = nullptr>
  std::size_t write(InputIt first_value, const InputIt &last_value, std::size_t offset = 0,
                    bool allow_dataset_resize = false, UnaryOperation op = identity());

  template <typename InputIt, typename UnaryOperation = identity,
            typename std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, InputIt>, int>
                * = nullptr>
  std::size_t append(InputIt first_value, const InputIt &last_value,
                     UnaryOperation op = identity());

  // Write single values
  template <typename N, typename = std::enable_if_t<std::is_arithmetic_v<N>>>
  std::size_t write(N buff, std::size_t offset = 0, bool allow_dataset_resize = false);
  std::size_t write(std::string buff, std::size_t offset = 0, bool allow_dataset_resize = false);
  std::size_t write(const GenericVariant &vbuff, std::size_t offset = 0,
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

  template <typename T>
  [[nodiscard]] auto begin(std::optional<std::ptrdiff_t> chunk_size_ = {}) const -> iterator<T>;
  template <typename T>
  [[nodiscard]] auto end(std::optional<std::ptrdiff_t> chunk_size_ = {}) const -> iterator<T>;

  template <typename T>
  [[nodiscard]] auto cbegin(std::optional<std::ptrdiff_t> chunk_size_ = {}) const -> iterator<T>;
  template <typename T>
  [[nodiscard]] auto cend(std::optional<std::ptrdiff_t> chunk_size_ = {}) const -> iterator<T>;

  template <typename T>
  [[nodiscard]] auto make_iterator_at_offset(std::size_t offset,
                                             std::optional<std::ptrdiff_t> chunk_size_ = {}) const
      -> iterator<T>;

  [[nodiscard]] static std::pair<std::string, std::string> parse_uri(std::string_view uri);

  template <typename T>
  [[nodiscard]] static auto lower_bound(iterator<T> first, iterator<T> last, const T &value,
                                        bool assume_uniform_distribution = false) -> iterator<T>;

 private:
  [[nodiscard]] const HighFive::Selection &select(std::size_t offset, std::size_t count = 1) const;
  [[nodiscard]] HighFive::Selection &select(std::size_t offset, std::size_t count = 1);

  [[nodiscard]] static HighFive::DataSet create_fixed_str_dataset(
      RootGroup &root_grp, std::string_view path, std::size_t max_str_length, std::size_t max_dim,
      const HighFive::DataSetAccessProps &aprops, const HighFive::DataSetCreateProps &cprops);

  [[noreturn]] void throw_out_of_range_excp(std::size_t offset) const;
  [[noreturn]] void throw_out_of_range_excp(std::size_t offset, std::size_t n) const;

  [[nodiscard]] HighFive::DataType get_h5type() const;
  [[nodiscard]] static std::size_t get_chunk_size(const HighFive::DataSet &dset);

  template <typename T>
  std::size_t read(T *buffer, std::size_t buff_size, std::size_t offset) const;

 public:
  template <typename T>
  class iterator {
    friend Dataset;
    mutable internal::COWChunk<T> _buffer{};
    std::shared_ptr<const Dataset> _dset{};
    std::uint32_t _chunk_size{};
    std::size_t _h5_offset{};
    std::size_t _h5_size{};

    explicit iterator(Dataset dset, std::optional<std::ptrdiff_t> chunk_size_ = {},
                      std::size_t h5_offset = 0, bool init = true);
    explicit iterator(std::shared_ptr<const Dataset> dset,
                      std::optional<std::ptrdiff_t> chunk_size_ = {}, std::size_t h5_offset = 0,
                      bool init = true);

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = T;
    using pointer = value_type *;
    using reference = value_type &;
    using iterator_category = std::random_access_iterator_tag;

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
    auto operator+=(difference_type i) -> iterator &;
    [[nodiscard]] auto operator+(difference_type i) const -> iterator;

    auto operator--() -> iterator &;
    auto operator--(int) -> iterator;
    auto operator-=(difference_type i) -> iterator &;
    [[nodiscard]] auto operator-(difference_type i) const -> iterator;
    [[nodiscard]] auto operator-(const iterator &other) const -> difference_type;

    template <typename I>
    auto seek(I offset) -> iterator &;
    [[nodiscard]] constexpr std::size_t h5_offset() const noexcept;
    [[nodiscard]] auto buffer() const -> const internal::COWChunk<T> &;

    [[nodiscard]] constexpr std::size_t chunk_size() const noexcept;
    [[nodiscard]] constexpr const Dataset &dataset() const noexcept;

   private:
    void read_chunk_at_offset(std::size_t new_offset, bool forward = true) const;
    [[nodiscard]] bool buffer_is_outdated() const noexcept;

    [[nodiscard]] static auto make_end_iterator(Dataset dset,
                                                std::optional<std::ptrdiff_t> chunk_size_ = {})
        -> iterator;
    [[nodiscard]] static auto make_end_iterator(std::shared_ptr<const Dataset> dset,
                                                std::optional<std::ptrdiff_t> chunk_size_ = {})
        -> iterator;
    [[nodiscard]] static std::uint32_t compute_chunk_size(
        const std::shared_ptr<const Dataset> &dset, std::optional<std::ptrdiff_t> chunk_size_);
    void bound_check(std::ptrdiff_t i = 0, bool close_interval = false) const noexcept;
  };
};
HICTK_DISABLE_WARNING_POP

using DatasetMap = phmap::flat_hash_map<std::string, Dataset>;

}  // namespace hictk::cooler

#include "./impl/dataset_accessors_impl.hpp"  // NOLINT
#include "./impl/dataset_impl.hpp"            // NOLINT
#include "./impl/dataset_iterator_impl.hpp"   // NOLINT
#include "./impl/dataset_read_impl.hpp"       // NOLINT
#include "./impl/dataset_write_impl.hpp"      // NOLINT
