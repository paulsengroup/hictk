// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
// clang-format off
#include "hictk/suppress_warnings.hpp"
// clang-format on
DISABLE_WARNING_PUSH
DISABLE_WARNING_NULL_DEREF
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
DISABLE_WARNING_POP

#include <initializer_list>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/numeric_variant.hpp"
#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

using DefaultPixelT = std::int32_t;

class File;
class MultiResFile;
class SingleCellFile;

struct Attributes {
  friend File;

  // Mandatory attributes
  std::uint32_t bin_size{0};
  std::optional<std::string> bin_type{"fixed"};  // Mandatory in v3
  std::string format{COOL_MAGIC};
  std::uint8_t format_version{3};
  std::optional<std::string> storage_mode{"symmetric-upper"};  // Mandatory in v3

  // Reserved attributes
  std::optional<std::string> creation_date{generate_creation_date()};
  std::optional<std::string> generated_by{HICTK_VERSION_STRING_LONG};
  std::optional<std::string> assembly{"unknown"};
  std::optional<std::string> metadata{"{}"};

  // Optional but common
  std::optional<std::string> format_url{"https://github.com/open2c/cooler"};
  std::optional<std::int64_t> nbins{0};
  std::optional<std::int32_t> nchroms{0};
  std::optional<std::int64_t> nnz{0};
  using SumVar = std::variant<double, std::int64_t>;
  std::optional<SumVar> sum{std::int64_t(0)};
  std::optional<SumVar> cis{std::int64_t(0)};

  template <typename PixelT = DefaultPixelT,
            typename = std::enable_if_t<std::is_arithmetic_v<PixelT>>>
  [[nodiscard]] static Attributes init(std::uint32_t bin_size_);
  [[nodiscard]] static Attributes init_empty() noexcept;
  [[nodiscard]] bool operator==(const Attributes &other) const noexcept;
  [[nodiscard]] bool operator!=(const Attributes &other) const noexcept;
  [[nodiscard]] static std::string generate_creation_date();

 private:
  // Use the init factory methods to construct an Attribute object
  Attributes() = default;
};

class File {
  friend MultiResFile;
  friend SingleCellFile;
  using NumericVariant = hictk::internal::NumericVariant;
  unsigned int _mode{HighFive::File::ReadOnly};
  RootGroup _root_group{};
  GroupMap _groups{};
  DatasetMap _datasets{};
  mutable balancing::WeightMap _weights{};
  mutable balancing::WeightMap _weights_scaled{};
  Attributes _attrs{Attributes::init(0)};
  NumericVariant _pixel_variant{};
  std::shared_ptr<const BinTable> _bins{};
  mutable std::shared_ptr<Index> _index{};
  bool _finalize{false};

  // Constructors are private. Cooler files are opened using factory methods
  File(RootGroup entrypoint, unsigned int mode, std::size_t cache_size_bytes, double w0,
       bool validate);

  template <typename PixelT>
  File(RootGroup entrypoint, Reference chroms, PixelT pixel, Attributes attributes,
       std::size_t cache_size_bytes, double w0);
  // Ctor for SingleCellCooler
  template <typename PixelT>
  File(RootGroup entrypoint, PixelT pixel, Attributes attributes, std::size_t cache_size_bytes,
       double w0);

 public:
  using QUERY_TYPE = hictk::GenomicInterval::Type;

  File() = default;
  File(const File &other) = delete;
  File(File &&other) noexcept(noexcept_move_ctor()) = default;  // NOLINT

  // Simple constructor. Open file in read-only mode. Automatically detects pixel count type
  explicit File(std::string_view uri, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
                bool validate = true);
  explicit File(RootGroup entrypoint, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
                bool validate = true);

  [[nodiscard]] static File open_random_access(
      std::string_view uri, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
      bool validate = true);
  [[nodiscard]] static File open_read_once(std::string_view uri,
                                           std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
                                           bool validate = true);
  template <typename PixelT = DefaultPixelT>
  [[nodiscard]] static File create(std::string_view uri, const Reference &chroms,
                                   std::uint32_t bin_size, bool overwrite_if_exists = false,
                                   Attributes attributes = Attributes::init<PixelT>(0),
                                   std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4);

  [[nodiscard]] static File open_random_access(
      RootGroup entrypoint, std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
      bool validate = true);
  [[nodiscard]] static File open_read_once(RootGroup entrypoint,
                                           std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE,
                                           bool validate = true);
  template <typename PixelT = DefaultPixelT>
  [[nodiscard]] static File create(RootGroup entrypoint, const Reference &chroms,
                                   std::uint32_t bin_size,
                                   Attributes attributes = Attributes::init<PixelT>(0),
                                   std::size_t cache_size_bytes = DEFAULT_HDF5_CACHE_SIZE * 4);

  ~File() noexcept;

  File &operator=(const File &other) = delete;
  File &operator=(File &&other) noexcept(noexcept_move_assignment_op()) = default;  // NOLINT

  [[nodiscard]] explicit operator bool() const noexcept;

  void close();

  [[nodiscard]] std::string uri() const;
  [[nodiscard]] std::string hdf5_path() const;
  [[nodiscard]] std::string path() const;

  [[nodiscard]] auto chromosomes() const noexcept -> const Reference &;
  [[nodiscard]] auto bins() const noexcept -> const BinTable &;
  [[nodiscard]] auto bins_ptr() const noexcept -> std::shared_ptr<const BinTable>;

  [[nodiscard]] std::uint32_t bin_size() const noexcept;
  [[nodiscard]] std::uint64_t nbins() const;
  [[nodiscard]] std::uint64_t nchroms() const;
  [[nodiscard]] std::uint64_t nnz() const;

  [[nodiscard]] auto attributes() const noexcept -> const Attributes &;
  [[nodiscard]] auto group(std::string_view group_name) -> Group &;
  [[nodiscard]] auto dataset(std::string_view dataset_name) -> Dataset &;
  [[nodiscard]] auto group(std::string_view group_name) const -> const Group &;
  [[nodiscard]] auto dataset(std::string_view dataset_name) const -> const Dataset &;

  [[nodiscard]] const NumericVariant &pixel_variant() const noexcept;
  template <typename T>
  [[nodiscard]] bool has_pixel_of_type() const noexcept;

  [[nodiscard]] bool has_signed_pixels() const noexcept;
  [[nodiscard]] bool has_unsigned_pixels() const noexcept;
  [[nodiscard]] bool has_integral_pixels() const noexcept;
  [[nodiscard]] bool has_float_pixels() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(PixelIt first_pixel, PixelIt last_pixel, bool validate = false);

  template <typename N>
  [[nodiscard]] typename PixelSelector::iterator<N> begin(
      std::string_view weight_name = "NONE") const;
  template <typename N>
  [[nodiscard]] typename PixelSelector::iterator<N> end(
      std::string_view weight_name = "NONE") const;

  template <typename N>
  [[nodiscard]] typename PixelSelector::iterator<N> cbegin(
      std::string_view weight_name = "NONE") const;
  template <typename N>
  [[nodiscard]] typename PixelSelector::iterator<N> cend(
      std::string_view weight_name = "NONE") const;

  [[nodiscard]] PixelSelector fetch(std::shared_ptr<const balancing::Weights> weights) const;
  [[nodiscard]] PixelSelector fetch(std::string_view range,
                                    std::shared_ptr<const balancing::Weights> weights,
                                    QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  [[nodiscard]] PixelSelector fetch(std::string_view chrom_name, std::uint32_t start,
                                    std::uint32_t end,
                                    std::shared_ptr<const balancing::Weights> weights) const;

  [[nodiscard]] PixelSelector fetch(std::string_view range1, std::string_view range2,
                                    std::shared_ptr<const balancing::Weights> weights,
                                    QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  [[nodiscard]] PixelSelector fetch(std::string_view chrom1_name, std::uint32_t start1,
                                    std::uint32_t end1, std::string_view chrom2_name,
                                    std::uint32_t start2, std::uint32_t end2,
                                    std::shared_ptr<const balancing::Weights> weights) const;

  [[nodiscard]] PixelSelector fetch(
      const balancing::Method &normalization = balancing::Method::NONE()) const;
  [[nodiscard]] PixelSelector fetch(
      std::string_view range, const balancing::Method &normalization = balancing::Method::NONE(),
      QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  [[nodiscard]] PixelSelector fetch(
      std::string_view chrom_name, std::uint32_t start, std::uint32_t end,
      const balancing::Method &normalization = balancing::Method::NONE()) const;

  [[nodiscard]] PixelSelector fetch(
      std::string_view range1, std::string_view range2,
      const balancing::Method &normalization = balancing::Method::NONE(),
      QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  [[nodiscard]] PixelSelector fetch(
      std::string_view chrom1_name, std::uint32_t start1, std::uint32_t end1,
      std::string_view chrom2_name, std::uint32_t start2, std::uint32_t end2,
      const balancing::Method &normalization = balancing::Method::NONE()) const;
  [[nodiscard]] PixelSelector fetch(
      std::uint64_t first_bin, std::uint64_t last_bin,
      std::shared_ptr<const balancing::Weights> weights = nullptr) const;
  [[nodiscard]] PixelSelector fetch(
      std::uint64_t first_bin1, std::uint64_t last_bin1, std::uint64_t first_bin2,
      std::uint64_t last_bin2, std::shared_ptr<const balancing::Weights> weights = nullptr) const;

  std::shared_ptr<const balancing::Weights> read_weights(std::string_view normalization,
                                                         bool rescale = false) const;
  std::shared_ptr<const balancing::Weights> read_weights(std::string_view normalization,
                                                         balancing::Weights::Type type,
                                                         bool rescale = false) const;

  [[nodiscard]] bool has_normalization(std::string_view normalization) const;
  [[nodiscard]] bool has_normalization(const balancing::Method &normalization) const;
  [[nodiscard]] std::vector<balancing::Method> avail_normalizations() const;
  std::shared_ptr<const balancing::Weights> read_weights(const balancing::Method &normalization,
                                                         bool rescale = false) const;
  std::shared_ptr<const balancing::Weights> read_weights(const balancing::Method &normalization,
                                                         balancing::Weights::Type type,
                                                         bool rescale = false) const;

  bool purge_weights(std::string_view name = "");

  void flush();

  template <typename It>
  static void write_weights(std::string_view uri, std::string_view name, It first_weight,
                            It last_weight, bool overwrite_if_exists = false,
                            bool divisive = false);
  template <typename It>
  void write_weights(std::string_view name, It first_weight, It last_weight,
                     bool overwrite_if_exists = false, bool divisive = false);

  void validate_bins(bool full = false) const;

 private:
  [[nodiscard]] auto index() const noexcept -> const Index &;
  [[nodiscard]] auto index() noexcept -> Index &;

  [[nodiscard]] static HighFive::File open_file(std::string_view uri, unsigned int mode,
                                                bool validate);

  [[nodiscard]] static auto open_or_create_root_group(HighFive::File f, std::string_view uri)
      -> RootGroup;

  // Open/read groups, datasets and attributes
  [[nodiscard]] static auto open_root_group(const HighFive::File &f, std::string_view uri)
      -> RootGroup;
  [[nodiscard]] static auto open_groups(const RootGroup &root_grp) -> GroupMap;
  [[nodiscard]] static auto open_datasets(const RootGroup &root_grp, std::size_t cache_size_bytes,
                                          double w0) -> DatasetMap;
  [[nodiscard]] static auto read_standard_attributes(const RootGroup &root_grp,
                                                     bool initialize_missing = false) -> Attributes;

  // Create/write groups, datasets and attributes
  [[nodiscard]] static auto create_root_group(HighFive::File &f, std::string_view uri,
                                              bool write_sentinel_attr = true) -> RootGroup;
  [[nodiscard]] static auto create_groups(RootGroup &root_grp) -> GroupMap;
  [[nodiscard]] static auto create_groups(RootGroup &root_grp, Group chroms_grp, Group bins_grp)
      -> GroupMap;
  template <typename PixelT>
  [[nodiscard]] static auto create_datasets(RootGroup &root_grp, const Reference &chroms,
                                            std::size_t cache_size_bytes, double w0) -> DatasetMap;
  static void write_standard_attributes(RootGroup &root_grp, const Attributes &attributes,
                                        bool skip_sentinel_attr = true);

  [[nodiscard]] static auto import_chroms(const Dataset &chrom_names, const Dataset &chrom_sizes,
                                          bool missing_ok) -> Reference;

  [[nodiscard]] static Index init_index(const Dataset &chrom_offset_dset,
                                        const Dataset &bin_offset_dset,
                                        std::shared_ptr<const BinTable> bin_table,
                                        std::uint64_t expected_nnz, bool missing_ok);
  void read_index_chunk(std::initializer_list<Chromosome> chroms) const;

  template <typename PixelIt>
  static void append_bins(Dataset &bin1_dset, Dataset &bin2_dset, PixelIt first_pixel,
                          PixelIt last_pixel);
  template <typename PixelIt, typename N>
  static void append_counts(Dataset &dset, const BinTable &bins, PixelIt first_pixel,
                            PixelIt last_pixel, N &sum, N &cis_sum);

  template <typename PixelIt>
  void validate_pixels_before_append(PixelIt first_pixel, PixelIt last_pixel) const;
  template <typename PixelIt>
  void validate_thin_pixels_before_append(PixelIt first_pixel, PixelIt last_pixel) const;

  [[nodiscard]] static NumericVariant detect_pixel_type(const RootGroup &root_grp,
                                                        std::string_view path = "pixels/count");
  void write_attributes(bool skip_sentinel_attr = true);
  void write_chromosomes();

  template <typename ChromIt, typename UnaryOperation = identity,
            typename = std::enable_if_t<is_unary_operation_on_iterator<UnaryOperation, ChromIt>>>
  static void write_chromosomes(Dataset &name_dset, Dataset &size_dset, ChromIt first_chrom,
                                ChromIt last_chrom, UnaryOperation op = identity());

  void write_bin_table();
  static void write_bin_table(Dataset &chrom_dset, Dataset &start_dset, Dataset &end_dset,
                              const BinTable &bin_table);
  template <typename PixelIt>
  void update_indexes(PixelIt first_pixel, PixelIt last_pixel);

  void write_indexes();
  static void write_indexes(Dataset &chrom_offset_dset, Dataset &bin_offset_dset, const Index &idx);

  void finalize();

  static void write_sentinel_attr(HighFive::Group grp);
  [[nodiscard]] static bool check_sentinel_attr(const HighFive::Group &grp);
  void write_sentinel_attr();
  [[nodiscard]] bool check_sentinel_attr();

  [[nodiscard]] Bin get_last_bin_written() const;

  template <typename N, bool cis = false>
  void update_pixel_sum(N partial_sum);

  template <typename PixelT>
  void validate_pixel_type() const noexcept;

  // IMPORTANT: the private fetch() methods interpret queries as open-open
  [[nodiscard]] PixelSelector fetch(PixelCoordinates coord,
                                    std::shared_ptr<const balancing::Weights> weights) const;
  [[nodiscard]] PixelSelector fetch(PixelCoordinates coord1, PixelCoordinates coord2,
                                    std::shared_ptr<const balancing::Weights> weights) const;
};

namespace internal {
template <typename N>
bool read_optional(const RootGroup &root_grp, std::string_view key, N &buff, bool missing_ok);
}

}  // namespace hictk::cooler

#include "./impl/file_accessors_impl.hpp"      // NOLINT
#include "./impl/file_impl.hpp"                // NOLINT
#include "./impl/file_read_impl.hpp"           // NOLINT
#include "./impl/file_standard_attr_impl.hpp"  // NOLINT
#include "./impl/file_validation_impl.hpp"     // NOLINT
#include "./impl/file_write_impl.hpp"          // NOLINT
