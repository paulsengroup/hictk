// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/string_utils.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/type_traits.hpp"

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

inline PixelSelector File::fetch(std::shared_ptr<const balancing::Weights> weights) const {
  if (!weights) {
    weights = normalization_ptr(balancing::Method::NONE());
  }

  // clang-format off
  return PixelSelector(
      _index,
      dataset("pixels/bin1_id"),
      dataset("pixels/bin2_id"),
      dataset("pixels/count"),
      std::move(weights));
  // clang-format on
}

inline PixelSelector File::fetch(std::string_view range,
                                 std::shared_ptr<const balancing::Weights> weights,
                                 QUERY_TYPE query_type) const {
  const auto gi = query_type == QUERY_TYPE::BED
                      ? GenomicInterval::parse_bed(chromosomes(), range)
                      : GenomicInterval::parse_ucsc(chromosomes(), std::string{range});

  return fetch(PixelCoordinates{bins().at(gi)}, std::move(weights));
}

inline PixelSelector File::fetch(std::string_view chrom_name, std::uint32_t start,
                                 std::uint32_t end,
                                 std::shared_ptr<const balancing::Weights> weights) const {
  assert(start < end);

  return fetch(PixelCoordinates{bins().at(chrom_name, start),
                                bins().at(chrom_name, end - (std::min)(end, 1U))},
               std::move(weights));
}

inline PixelSelector File::fetch(PixelCoordinates coord,
                                 std::shared_ptr<const balancing::Weights> weights) const {
  const auto &current_chrom = coord.bin1.chrom();
  const auto &next_chrom = chromosomes().at(
      std::min(static_cast<std::uint32_t>(chromosomes().size() - 1), coord.bin1.chrom().id() + 1));
  read_index_chunk({current_chrom, next_chrom});

  if (!weights) {
    weights = normalization_ptr(balancing::Method::NONE());
  }

  // clang-format off
  return PixelSelector(_index,
                       dataset("pixels/bin1_id"),
                       dataset("pixels/bin2_id"),
                       dataset("pixels/count"),
                       std::move(coord),
                       std::move(weights)
  );
  // clang-format on
}

inline PixelSelector File::fetch(std::string_view range1, std::string_view range2,
                                 std::shared_ptr<const balancing::Weights> weights,
                                 QUERY_TYPE query_type) const {
  if (range1 == range2) {
    return fetch(range1, std::move(weights), query_type);
  }

  const auto gi1 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(chromosomes(), range1)
                       : GenomicInterval::parse_ucsc(chromosomes(), std::string{range1});

  const auto gi2 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(chromosomes(), range2)
                       : GenomicInterval::parse_ucsc(chromosomes(), std::string{range2});

  return fetch(PixelCoordinates{bins().at(gi1)}, PixelCoordinates{bins().at(gi2)},
               std::move(weights));
}

inline PixelSelector File::fetch(std::string_view chrom1, std::uint32_t start1, std::uint32_t end1,
                                 std::string_view chrom2, std::uint32_t start2, std::uint32_t end2,
                                 std::shared_ptr<const balancing::Weights> weights) const {
  assert(start1 < end1);
  assert(start2 < end2);

  PixelCoordinates coord1{bins().at(chrom1, start1),
                          bins().at(chrom1, end1 - (std::min)(end1, 1U))};
  PixelCoordinates coord2{bins().at(chrom2, start2),
                          bins().at(chrom2, end2 - (std::min)(end2, 1U))};

  return fetch(coord1, coord2, std::move(weights));
}
inline PixelSelector File::fetch(const balancing::Method &normalization_) const {
  return fetch(normalization_ptr(normalization_));
}
inline PixelSelector File::fetch(std::string_view range, const balancing::Method &normalization_,
                                 QUERY_TYPE query_type) const {
  return fetch(range, normalization_ptr(normalization_), query_type);
}
inline PixelSelector File::fetch(std::string_view chrom_name, std::uint32_t start,
                                 std::uint32_t end, const balancing::Method &normalization_) const {
  return fetch(chrom_name, start, end, normalization_ptr(normalization_));
}

inline PixelSelector File::fetch(std::string_view range1, std::string_view range2,
                                 const balancing::Method &normalization_,
                                 QUERY_TYPE query_type) const {
  return fetch(range1, range2, normalization_ptr(normalization_), query_type);
}
inline PixelSelector File::fetch(std::string_view chrom1_name, std::uint32_t start1,
                                 std::uint32_t end1, std::string_view chrom2_name,
                                 std::uint32_t start2, std::uint32_t end2,
                                 const balancing::Method &normalization_) const {
  return fetch(chrom1_name, start1, end1, chrom2_name, start2, end2,
               normalization_ptr(normalization_));
}

inline PixelSelector File::fetch(std::uint64_t first_bin, std::uint64_t last_bin,
                                 std::shared_ptr<const balancing::Weights> weights) const {
  return fetch(first_bin, last_bin, first_bin, last_bin, std::move(weights));
}

inline PixelSelector File::fetch(std::uint64_t first_bin1, std::uint64_t last_bin1,
                                 std::uint64_t first_bin2, std::uint64_t last_bin2,
                                 std::shared_ptr<const balancing::Weights> weights) const {
  PixelCoordinates coord1{bins().at(first_bin1), bins().at(last_bin1)};
  PixelCoordinates coord2{bins().at(first_bin2), bins().at(last_bin2)};

  return fetch(coord1, coord2, std::move(weights));
}

inline PixelSelector File::fetch(PixelCoordinates coord1, PixelCoordinates coord2,
                                 std::shared_ptr<const balancing::Weights> weights) const {
  const auto &current_chrom = coord1.bin1.chrom();
  const auto &next_chrom = chromosomes().at(
      std::min(static_cast<std::uint32_t>(chromosomes().size() - 1), coord1.bin1.chrom().id() + 1));
  read_index_chunk({current_chrom, next_chrom});

  if (!weights) {
    weights = normalization_ptr(balancing::Method::NONE());
  }

  // clang-format off
  return PixelSelector(_index,
                       dataset("pixels/bin1_id"),
                       dataset("pixels/bin2_id"),
                       dataset("pixels/count"),
                       std::move(coord1),
                       std::move(coord2),
                       std::move(weights)
  );
  // clang-format on
}

inline bool File::has_normalization(std::string_view normalization_) const {
  return has_normalization(balancing::Method{normalization_});
}
inline const balancing::Weights &File::normalization(std::string_view normalization_,
                                                     bool rescale) const {
  return *normalization_ptr(balancing::Method{normalization_}, rescale);
}
inline const balancing::Weights &File::normalization(std::string_view normalization_,
                                                     balancing::Weights::Type type,
                                                     bool rescale) const {
  return *normalization_ptr(balancing::Method{normalization_}, type, rescale);
}

inline const balancing::Weights &File::normalization(const balancing::Method &normalization_,
                                                     bool rescale) const {
  return *normalization_ptr(normalization_, rescale);
}
inline std::shared_ptr<const balancing::Weights> File::normalization_ptr(
    std::string_view normalization_, bool rescale) const {
  return normalization_ptr(balancing::Method{normalization_}, rescale);
}
inline std::shared_ptr<const balancing::Weights> File::normalization_ptr(
    std::string_view normalization_, balancing::Weights::Type type, bool rescale) const {
  return normalization_ptr(balancing::Method{normalization_}, type, rescale);
}

inline std::shared_ptr<const balancing::Weights> File::normalization_ptr(
    const balancing::Method &normalization_, bool rescale) const {
  return normalization_ptr(normalization_, balancing::Weights::Type::INFER, rescale);
}

inline std::shared_ptr<const balancing::Weights> File::normalization_ptr(
    const balancing::Method &normalization_, balancing::Weights::Type type, bool rescale) const {
  if (!rescale) {
    if (const auto it = _weights.find(normalization_.to_string()); it != _weights.end()) {
      return it->second;
    }
  } else {
    if (const auto it = _weights_scaled.find(normalization_.to_string());
        it != _weights_scaled.end()) {
      return it->second;
    }
  }

  if (normalization_ == "NONE") {
    auto weights = std::make_shared<const balancing::Weights>(
        1.0, bins().size(), balancing::Weights::Type::MULTIPLICATIVE);
    if (rescale) {
      _weights_scaled.emplace(normalization_.to_string(), weights);
    } else {
      _weights.emplace(normalization_.to_string(), weights);
    }

    return weights;
  }

  const auto dset_path = fmt::format(FMT_STRING("{}/{}"), _groups.at("bins").group.getPath(),
                                     normalization_.to_string());

  if (!_root_group().exist(dset_path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to read \"{}\" weights: dataset \"{}\" does not exist"),
                    normalization_.to_string(), dset_path));
  }

  Dataset dset{
      _root_group, dset_path,
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

  auto weights = std::make_shared<balancing::Weights>(dset.read_all<std::vector<double>>(), type);
  if (!rescale) {
    _weights.emplace(normalization_.to_string(), weights);
    return weights;
  }

  if (!dset.has_attribute("scale")) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to read scaling factors from {}"), dset.hdf5_path()));
  }

  const auto cis_only =
      dset.has_attribute("cis_only") ? dset.read_attribute<bool>("cis_only") : false;

  if (cis_only) {
    std::vector<double> scaling_factors;
    dset.read_attribute("scale", scaling_factors);

    const auto bin_offsets = bins().num_bin_prefix_sum();

    assert(!bin_offsets.empty());
    if (bin_offsets.size() - 1 != scaling_factors.size()) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("failed to read weights from \"{}\": expected {} scale value(s), found {}"),
          dset.uri(), bin_offsets.size() - 1, scaling_factors.size()));
    }

    weights->rescale(scaling_factors, bin_offsets);
  } else {
    weights->rescale(dset.read_attribute<double>("scale"));
  }

  _weights_scaled.emplace(normalization_.to_string(), weights);
  return weights;
}

inline bool File::purge_weights(std::string_view name) {
  if (_weights.empty()) {
    return false;
  }
  if (name == "") {
    _weights.clear();
    return true;
  }
  return _weights.erase(std::string{name});
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

inline auto File::open_datasets(const RootGroup &root_grp, std::size_t cache_size_bytes,
                                double w0) -> DatasetMap {
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
inline bool read_optional(const RootGroup &root_grp, std::string_view key, N &buff,
                          bool missing_ok) {
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

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_UNREACHABLE_CODE
inline auto File::read_standard_attributes(const RootGroup &root_grp,
                                           bool initialize_missing) -> Attributes {
  auto attrs = initialize_missing ? Attributes::init(0) : Attributes::init_empty();
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
  read_or_throw("format", attrs.format);

  // Read mandatory attributes for Cooler v3
  auto missing_ok = attrs.format_version < 3;

  std::optional<std::string> bin_type{"fixed"};
  internal::read_optional(root_grp, "bin-type", bin_type, missing_ok);
  attrs.bin_type = bin_type.value() == "fixed" ? BinTable::Type::fixed : BinTable::Type::variable;
  if (attrs.bin_type == BinTable::Type::fixed) {
    read_or_throw("bin-size", attrs.bin_size);
  }
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
HICTK_DISABLE_WARNING_POP

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

inline BinTable File::init_bin_table(const DatasetMap &dsets, BinTable::Type bin_type,
                                     std::uint32_t bin_size) {
  auto chroms = import_chroms(dsets.at("chroms/name"), dsets.at("chroms/length"), false);
  if (bin_type == BinTable::Type::fixed) {
    return {std::move(chroms), bin_size};
  }
  assert(bin_type == BinTable::Type::variable);
  assert(bin_size == 0);

  return {std::move(chroms), dsets.at("bins/start").read_all<std::vector<std::uint32_t>>(),
          dsets.at("bins/end").read_all<std::vector<std::uint32_t>>()};
}

inline Index File::init_index(const Dataset &chrom_offset_dset, const Dataset &bin_offset_dset,
                              std::shared_ptr<const BinTable> bin_table, std::uint64_t expected_nnz,
                              bool missing_ok) {
  assert(bin_table);
  try {
    if (bin_offset_dset.empty()) {
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

    auto first = bin_offset_dset.begin<std::uint64_t>();
    auto last = bin_offset_dset.end<std::uint64_t>();

    assert(first != last);
    if (const auto offset = *first; offset != 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("{} is corrupted: first offset should be 0, found {}"),
                      bin_offset_dset.hdf5_path(), offset));
    }
    if (const auto offset = *(last - 1); offset != expected_nnz) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("{} is corrupted: last offset should be {}, found {}"),
                      bin_offset_dset.hdf5_path(), expected_nnz, offset));
    }
    std::vector<std::uint64_t> chrom_offsets{};
    for (const auto &chrom_offset :
         internal::import_chrom_offsets(chrom_offset_dset, bin_table->chromosomes().size() + 1)) {
      const auto bin1_offset = bin_offset_dset.read<std::uint64_t>(chrom_offset);
      chrom_offsets.push_back(bin1_offset);
    }

    return Index{bin_table, chrom_offsets, expected_nnz, false};
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to initialize index for cooler at URI: \"{}\": {}"),
                    bin_offset_dset.get_parent().uri(), e.what()));
  }
}

inline void File::read_index_chunk(std::initializer_list<Chromosome> chroms) const {
  assert(_index);
  try {
    for (const auto &chrom : chroms) {
      if (_index->size(chrom.id()) != 1) {
        continue;
      }

      auto chrom_offset_dset = dataset("indexes/chrom_offset");
      auto bin_offset_dset = dataset("indexes/bin1_offset");
      const auto chrom_offsets =
          internal::import_chrom_offsets(chrom_offset_dset, chromosomes().size() + 1);

      auto offset1 = static_cast<std::ptrdiff_t>(chrom_offsets[chrom.id()]);
      auto offset2 = static_cast<std::ptrdiff_t>(chrom_offsets[chrom.id() + 1]);
      auto first = bin_offset_dset.begin<std::uint64_t>() + offset1;
      auto last = bin_offset_dset.begin<std::uint64_t>() + offset2;
      _index->set(chrom, {first, last});

      try {
        _index->validate(chrom);
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(FMT_STRING("index validation failed: {}"), e.what()));
      }
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to import indexes for cooler at URI: \"{}\": {}"),
                    dataset("indexes/bin1_offset").get_parent().uri(), e.what()));
  }
}

inline bool File::check_sentinel_attr() { return File::check_sentinel_attr(_root_group()); }

inline Bin File::get_last_bin_written() const {
  const auto &dset = dataset("pixels/bin1_id");
  if (dset.empty()) {
    return bins().at(0);
  }
  const auto bin1_id = dset.read_last<std::uint64_t>();
  return bins().at(bin1_id);
}

}  // namespace hictk::cooler
