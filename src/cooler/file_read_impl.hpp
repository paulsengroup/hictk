// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/string_utils.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/type_pretty_printer.hpp"

namespace hictk::cooler {

inline bool File::check_sentinel_attr(const HighFive::Group &grp) {
  const auto generated_by_v = Attribute::read(grp, "generated-by", true);
  if (const auto *generated_by = std::get_if<std::string>(&generated_by_v);
      !generated_by || generated_by->find("hictk") == std::string::npos) {
    return false;
  }

  using T = remove_cvref_t<decltype(internal::SENTINEL_ATTR_VALUE)>;
  const auto sentinel_v = Attribute::read(grp, internal::SENTINEL_ATTR_NAME, true);
  const auto *sentinel = std::get_if<T>(&sentinel_v);

  return static_cast<bool>(sentinel) && *sentinel == internal::SENTINEL_ATTR_VALUE;
}

namespace internal {
[[nodiscard]] inline std::vector<std::uint64_t> import_chrom_offsets(const Dataset &dset,
                                                                     std::size_t expected_size) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto offsets = dset.read_all<std::vector<std::uint64_t>>();
  try {
    if (offsets.size() != expected_size) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("expected {} offsets, found {}"), expected_size, offsets.size()));
    }
    if (offsets.front() != 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("first offset should be 0, found {}"), offsets.front()));
    }
    if (!std::is_sorted(offsets.begin(), offsets.end())) {
      throw std::runtime_error("offsets are not in ascending order");
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("failed to import offsets from {}: {}"), dset.uri(), e.what()));
  }

  return offsets;
}
}  // namespace internal

template <std::size_t CHUNK_SIZE>
inline PixelSelector<CHUNK_SIZE> File::fetch(
    std::shared_ptr<const balancing::Weights> weights) const {
  // clang-format off
  return PixelSelector<CHUNK_SIZE>(
      this->_index,
      this->dataset("pixels/bin1_id"),
      this->dataset("pixels/bin2_id"),
      this->dataset("pixels/count"),
      weights);
  // clang-format on
}

template <std::size_t CHUNK_SIZE>
inline PixelSelector<CHUNK_SIZE> File::fetch(std::string_view query,
                                             std::shared_ptr<const balancing::Weights> weights,
                                             QUERY_TYPE query_type) const {
  const auto gi = query_type == QUERY_TYPE::BED
                      ? GenomicInterval::parse_bed(this->chromosomes(), query)
                      : GenomicInterval::parse_ucsc(this->chromosomes(), std::string{query});

  return this->fetch<CHUNK_SIZE>(PixelCoordinates{this->bins().at(gi)}, std::move(weights));
}

template <std::size_t CHUNK_SIZE>
inline PixelSelector<CHUNK_SIZE> File::fetch(
    std::string_view chrom_name, std::uint32_t start, std::uint32_t end,
    std::shared_ptr<const balancing::Weights> weights) const {
  assert(start < end);

  return this->fetch<CHUNK_SIZE>(
      PixelCoordinates{this->bins().at(chrom_name, start),
                       this->bins().at(chrom_name, end - (std::min)(end, 1U))},
      std::move(weights));
}

template <std::size_t CHUNK_SIZE>
inline PixelSelector<CHUNK_SIZE> File::fetch(
    PixelCoordinates coord, std::shared_ptr<const balancing::Weights> weights) const {
  // clang-format off
  return PixelSelector<CHUNK_SIZE>(this->_index,
                                   this->dataset("pixels/bin1_id"),
                                   this->dataset("pixels/bin2_id"),
                                   this->dataset("pixels/count"),
                                   std::move(coord),
                                   std::move(weights)
  );
  // clang-format on
}

template <std::size_t CHUNK_SIZE>
inline PixelSelector<CHUNK_SIZE> File::fetch(std::string_view range1, std::string_view range2,
                                             std::shared_ptr<const balancing::Weights> weights,
                                             QUERY_TYPE query_type) const {
  if (range1 == range2) {
    return this->fetch<CHUNK_SIZE>(range1);
  }

  const auto gi1 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(this->chromosomes(), range1)
                       : GenomicInterval::parse_ucsc(this->chromosomes(), std::string{range1});

  const auto gi2 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(this->chromosomes(), range2)
                       : GenomicInterval::parse_ucsc(this->chromosomes(), std::string{range2});

  return this->fetch<CHUNK_SIZE>(PixelCoordinates{this->bins().at(gi1)},
                                 PixelCoordinates{this->bins().at(gi2)}, std::move(weights));
}

template <std::size_t CHUNK_SIZE>
inline PixelSelector<CHUNK_SIZE> File::fetch(
    std::string_view chrom1, std::uint32_t start1, std::uint32_t end1, std::string_view chrom2,
    std::uint32_t start2, std::uint32_t end2,
    std::shared_ptr<const balancing::Weights> weights) const {
  assert(start1 < end1);
  assert(start2 < end2);
  // clang-format off
  return PixelSelector<CHUNK_SIZE>(this->_index,
                                   this->dataset("pixels/bin1_id"),
                                   this->dataset("pixels/bin2_id"),
                                   this->dataset("pixels/count"),
                                   PixelCoordinates{this->bins().at(chrom1, start1),
                                                    this->bins().at(chrom1, end1 - (std::min)(end1, 1U))},
                                   PixelCoordinates{this->bins().at(chrom2, start2),
                                                    this->bins().at(chrom2, end2 - (std::min)(end2, 1U))},
                                   std::move(weights)
  );
  // clang-format on
}

template <std::size_t CHUNK_SIZE>
inline PixelSelector<CHUNK_SIZE> File::fetch(
    PixelCoordinates coord1, PixelCoordinates coord2,
    std::shared_ptr<const balancing::Weights> weights) const {
  // clang-format off
  return PixelSelector<CHUNK_SIZE>(this->_index,
                                   this->dataset("pixels/bin1_id"),
                                   this->dataset("pixels/bin2_id"),
                                   this->dataset("pixels/count"),
                                   std::move(coord1),
                                   std::move(coord2),
                                   std::move(weights)
  );
  // clang-format on
}

inline bool File::has_weights(std::string_view name) const {
  const auto dset_path =
      fmt::format(FMT_STRING("{}/{}"), this->_groups.at("bins").group.getPath(), name);
  if (this->_weights.contains(dset_path)) {
    return true;
  }

  return this->_root_group().exist(dset_path);
}

inline std::shared_ptr<const balancing::Weights> File::read_weights(std::string_view name) const {
  if (name == "NONE") {
    return nullptr;
  }
  if (name.empty()) {
    throw std::runtime_error("weight dataset name is empty");
  }

  return this->read_weights(name, balancing::Weights::infer_type(name));
}

inline std::shared_ptr<const balancing::Weights> File::read_weights(
    std::string_view name, balancing::Weights::Type type) const {
  if (name.empty()) {
    throw std::runtime_error("weight dataset name is empty");
  }


  const auto dset_path =
      fmt::format(FMT_STRING("{}/{}"), this->_groups.at("bins").group.getPath(), name);
  if (const auto it = _weights.find(dset_path); it != _weights.end()) {
    return it->second;
  }

  if (!this->_root_group().exist(dset_path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to read \"{}\" weights: dataset \"{}\" does not exist"),
                    name, dset_path));
  }

  Dataset dset{
      this->_root_group, dset_path,
      Dataset::init_access_props(DEFAULT_HDF5_CHUNK_SIZE, DEFAULT_HDF5_DATASET_CACHE_SIZE, 1.0)};

  if (type == balancing::Weights::Type::INFER || type == balancing::Weights::Type::UNKNOWN) {
    if (dset.has_attribute("divisive_weights")) {
      type = dset.read_attribute<bool>("divisive_weights")
                 ? balancing::Weights::Type::DIVISIVE
                 : balancing::Weights::Type::MULTIPLICATIVE;
    } else {
      type = balancing::Weights::infer_type(dset.name());
      if (type == balancing::Weights::Type::UNKNOWN) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("unable to infer type for \"{}\" weights"), dset.uri()));
      }
    }
  }

  const auto node = this->_weights.emplace(
      name, std::make_shared<const balancing::Weights>(dset.read_all<std::vector<double>>(), type));
  return node.first->second;
}

inline bool File::purge_weights(std::string_view name) {
  if (this->_weights.empty()) {
    return false;
  }
  if (name == "") {
    this->_weights.clear();
    return true;
  }
  return this->_weights.erase(std::string{name});
}

inline auto File::open_root_group(const HighFive::File &f, std::string_view uri) -> RootGroup {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  RootGroup grp{f.getGroup(parse_cooler_uri(uri).group_path)};
  if (File::check_sentinel_attr(grp())) {
    throw std::runtime_error("file was not properly closed");
  }
  return grp;
}

inline auto File::open_groups(const RootGroup &root_grp) -> GroupMap {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  GroupMap groups(MANDATORY_GROUP_NAMES.size() + 1);
  groups.emplace(root_grp.hdf5_path(), Group{root_grp, root_grp()});

  std::transform(MANDATORY_GROUP_NAMES.begin(), MANDATORY_GROUP_NAMES.end(),
                 std::inserter(groups, groups.begin()), [&root_grp](const auto group_name) {
                   const auto name = std::string{group_name};
                   auto group_obj = root_grp().getGroup(std::string{group_name});

                   return std::make_pair(name, Group{root_grp, group_obj});
                 });
  return groups;
}

inline auto File::open_datasets(const RootGroup &root_grp, std::size_t cache_size_bytes, double w0)
    -> DatasetMap {
  DatasetMap datasets(MANDATORY_DATASET_NAMES.size());

  const std::size_t num_pixel_datasets = 3;
  const std::size_t num_read_once_dataset = MANDATORY_DATASET_NAMES.size() - num_pixel_datasets;

  const std::size_t read_once_cache_size = DEFAULT_HDF5_DATASET_CACHE_SIZE;
  const std::size_t pixel_dataset_cache_size =
      (cache_size_bytes - (read_once_cache_size * num_read_once_dataset)) / num_pixel_datasets;

  const auto default_aprop =
      Dataset::init_access_props(DEFAULT_HDF5_CHUNK_SIZE, read_once_cache_size, 1.0);
  const auto pixels_aprop = Dataset::init_access_props(
      DEFAULT_HDF5_CHUNK_SIZE, ((std::max)(read_once_cache_size, pixel_dataset_cache_size)), w0);

  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto open_dataset = [&](const auto dataset_uri) {
    return std::make_pair(
        std::string{dataset_uri},
        Dataset{
            root_grp, dataset_uri,
            hictk::internal::starts_with(dataset_uri, "pixels") ? pixels_aprop : default_aprop});
  };

  std::transform(MANDATORY_DATASET_NAMES.begin(), MANDATORY_DATASET_NAMES.end(),
                 std::inserter(datasets, datasets.begin()), open_dataset);

  return datasets;
}

namespace internal {
template <typename N>
bool read_optional(const RootGroup &root_grp, std::string_view key, N &buff, bool missing_ok) {
  if (!Attribute::exists(root_grp(), key) && missing_ok) {
    return false;
  }

  try {
    using T = remove_cvref_t<decltype(*buff)>;
    buff = Attribute::read<T>(root_grp(), key);
    return true;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"), key,
                    root_grp().getPath(), e.what()));
  }
}

template <typename N>
bool read_sum_optional(const RootGroup &root_grp, std::string_view key, N &buff, bool missing_ok) {
  if (!Attribute::exists(root_grp(), key) && missing_ok) {
    return false;
  }

  try {
    auto sumv = Attribute::read(root_grp(), key);
    const auto ok = std::visit(
        [&](auto sum) {
          using T = remove_cvref_t<decltype(sum)>;
          if constexpr (std::is_integral_v<T>) {
            buff = conditional_static_cast<std::int64_t>(sum);
            return true;
          }
          if constexpr (std::is_floating_point_v<T>) {
            buff = conditional_static_cast<double>(sum);
            return true;
          }
          return false;
        },
        sumv);

    if (!ok) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("attribute \"{}{}\" does not have a numeric type"),
                      root_grp().getPath(), key));
    }

    return ok;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"), key,
                    root_grp().getPath(), e.what()));
  }
}

}  // namespace internal

DISABLE_WARNING_PUSH
DISABLE_WARNING_UNREACHABLE_CODE
inline auto File::read_standard_attributes(const RootGroup &root_grp, bool initialize_missing)
    -> StandardAttributes {
  auto attrs = initialize_missing ? StandardAttributes::init(0) : StandardAttributes::init_empty();
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT

  auto read_or_throw = [&](const auto &key, auto &buff) {
    using T = remove_cvref_t<decltype(buff)>;
    try {
      buff = Attribute::read<T>(root_grp(), key);
    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"),
                      key, root_grp().getPath(), e.what()));
    }
  };

  // Read mandatory attributes
  // We read format-version first because some attributes are mandatory only for cooler v3
  read_or_throw("format-version", attrs.format_version);
  read_or_throw("bin-size", attrs.bin_size);
  read_or_throw("format", attrs.format);

  // Read mandatory attributes for Cooler v3
  auto missing_ok = attrs.format_version < 3;
  internal::read_optional(root_grp, "bin-type", attrs.bin_type, missing_ok);
  internal::read_optional(root_grp, "storage-mode", attrs.storage_mode, missing_ok);

  // Try to read reserved attributes
  missing_ok = true;
  internal::read_optional(root_grp, "creation-date", attrs.creation_date, missing_ok);
  internal::read_optional(root_grp, "format-url", attrs.format_url, missing_ok);
  internal::read_optional(root_grp, "generated-by", attrs.generated_by, missing_ok);

  if (!internal::read_optional(root_grp, "genome-assembly", attrs.assembly, missing_ok)) {
    internal::read_optional(root_grp, "assembly", attrs.assembly, missing_ok);
  }

  internal::read_optional(root_grp, "metadata", attrs.metadata, missing_ok);

  // Try to read other common attributes
  internal::read_optional(root_grp, "nbins", attrs.nbins, missing_ok);
  internal::read_optional(root_grp, "nchroms", attrs.nchroms, missing_ok);
  internal::read_optional(root_grp, "nnz", attrs.nnz, missing_ok);

  internal::read_sum_optional(root_grp, "sum", attrs.sum, missing_ok);
  internal::read_sum_optional(root_grp, "cis", attrs.cis, missing_ok);

  return attrs;
}
DISABLE_WARNING_POP

inline auto File::import_chroms(const Dataset &chrom_names, const Dataset &chrom_sizes,
                                bool missing_ok) -> Reference {
  try {
    [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
    std::vector<std::string> names;
    std::vector<std::uint32_t> sizes;
    chrom_names.read_all(names);
    chrom_sizes.read_all(sizes);

    if (names.size() != sizes.size()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Cooler file \"{}\" appears to be corrupted: {} and "
                                 "{} shape mismatch: found {} name(s) and {} length(s)"),
                      chrom_names.file_name(), chrom_names.hdf5_path(), chrom_sizes.hdf5_path(),
                      names.size(), sizes.size()));
    }

    return Reference{names.begin(), names.end(), sizes.begin()};

  } catch ([[maybe_unused]] const HighFive::Exception &e) {
    if (missing_ok) {
      return {};
    }
    throw;
  }
}

inline Index File::import_indexes(const Dataset &chrom_offset_dset, const Dataset &bin_offset_dset,
                                  const Reference &chroms,
                                  std::shared_ptr<const BinTable> bin_table,
                                  std::uint64_t expected_nnz, bool missing_ok) {
  assert(bin_table);
  try {
    if (bin_offset_dset.empty()) {
      assert(chrom_offset_dset.empty());
      if (missing_ok) {
        return Index{bin_table};
      }
      throw std::runtime_error("index datasets are empty");
    }

    if (bin_offset_dset.size() != bin_table->size() + 1) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("failed to import offsets from {}: expected {} offsets, found {}"),
                      bin_offset_dset.hdf5_path(), bin_table->size() + 1, bin_offset_dset.size()));
    }

    const auto chrom_offsets = internal::import_chrom_offsets(chrom_offset_dset, chroms.size() + 1);

    Index idx{bin_table, expected_nnz};

    std::size_t bin_id = 0;
    std::for_each(bin_offset_dset.begin<std::uint64_t>(), bin_offset_dset.end<std::uint64_t>(),
                  [&](std::uint64_t offset) {
                    if (bin_id < bin_table->size()) {
                      idx.set_offset_by_bin_id(bin_id++, offset);
                    } else {
                      // Last bin
                      assert(bin_id == bin_table->size());
                    }
                  });

    try {
      idx.validate();
    } catch (const std::exception &e) {
      throw std::runtime_error(fmt::format(FMT_STRING("index validation failed: {}"), e.what()));
    }

    return idx;

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to import indexes for cooler at URI: \"{}\": {}"),
                    bin_offset_dset.get_parent().uri(), e.what()));
  }
}

inline bool File::check_sentinel_attr() { return File::check_sentinel_attr(this->_root_group()); }

inline Bin File::get_last_bin_written() const {
  const auto &dset = this->dataset("pixels/bin1_id");
  if (dset.empty()) {
    return this->bins().at(0);
  }
  const auto bin1_id = dset.read_last<std::uint64_t>();
  return this->bins().at(bin1_id);
}

}  // namespace hictk::cooler
