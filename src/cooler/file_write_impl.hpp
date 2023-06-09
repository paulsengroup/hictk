// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/uri.hpp"

namespace hictk {

template <typename PixelIt, typename>
inline void File::append_pixels(PixelIt first_pixel, PixelIt last_pixel, bool validate) {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);

  if constexpr (ndebug_not_defined()) {
    this->validate_pixel_type<T>();
  }

  this->update_indexes(first_pixel, last_pixel);

  if (validate) {
    this->validate_pixels_before_append(first_pixel, last_pixel);
  }

  this->dataset("pixels/bin1_id").append(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
    return pixel.coords.bin1.id();
  });

  this->dataset("pixels/bin2_id").append(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
    return pixel.coords.bin2.id();
  });

  T sum = 0;
  T cis_sum = 0;
  this->dataset("pixels/count").append(first_pixel, last_pixel, [&](const Pixel<T> &pixel) {
    if (pixel.count == 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Found pixel with 0 interactions: {}"), pixel.coords));
    }
    sum += pixel.count;
    if (pixel.coords.bin1.chrom().id() == pixel.coords.bin2.chrom().id()) {
      cis_sum += pixel.count;
    }
    return pixel.count;
  });

  this->_attrs.nnz = this->dataset("pixels/bin1_id").size();

  this->update_pixel_sum(sum);
  this->update_pixel_sum<T, true>(cis_sum);
}

inline void File::flush() { this->_fp->flush(); }

template <typename It>
inline void File::write_weights(std::string_view uri, std::string_view name, It first_weight,
                                It last_weight, bool overwrite_if_exists, bool divisive) {
  File(uri, HighFive::File::ReadWrite)
      .write_weights(name, first_weight, last_weight, overwrite_if_exists, divisive);
}

template <typename It>
inline void File::write_weights(std::string_view name, It first_weight, It last_weight,
                                bool overwrite_if_exists, bool divisive) {
  if (name.empty()) {
    throw std::runtime_error("weight name is empty");
  }

  if (this->_mode == HighFive::File::ReadOnly) {
    throw std::runtime_error("File::write_weights() was called on a file open in read-only mode");
  }

  const auto num_weights = std::distance(first_weight, last_weight);
  const auto expected_num_weights = static_cast<std::ptrdiff_t>(this->bins().size());
  if (num_weights != expected_num_weights) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Invalid weight shape, expected {} values, found {}"),
                    expected_num_weights, num_weights));
  }

  auto dset = [&, name_ = std::string{name}]() {
    // Return existing dataset
    auto &grp = this->group("bins").group;
    if (overwrite_if_exists && grp.exist(name_)) {
      return Dataset(this->_root_group, grp.getDataSet(name_));
    }

    // Create new dataset or throw
    const auto path = fmt::format(FMT_STRING("bins/{}"), name);
    typename std::iterator_traits<It>::value_type buff{};
    return Dataset(this->_root_group, path, buff, HighFive::DataSpace::UNLIMITED);
  }();

  dset.resize(static_cast<std::size_t>(std::distance(first_weight, last_weight)));
  dset.write(first_weight, last_weight);
  dset.write_attribute("divisive_weights", std::uint8_t(divisive), overwrite_if_exists);
}

inline auto File::create_root_group(HighFive::File &f, std::string_view uri,
                                    bool write_sentinel_attr) -> RootGroup {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto grp = f.createGroup(parse_cooler_uri(uri).group_path);
  if (write_sentinel_attr) {
    Attribute::write(grp, internal::SENTINEL_ATTR_NAME, internal::SENTINEL_ATTR_VALUE);
    f.flush();
  }

  return {grp};
}

inline auto File::create_groups(RootGroup &root_grp) -> GroupMap {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  GroupMap groups(MANDATORY_GROUP_NAMES.size() + 1);
  groups.emplace(root_grp.hdf5_path(), Group{root_grp, root_grp()});

  std::transform(MANDATORY_GROUP_NAMES.begin(), MANDATORY_GROUP_NAMES.end(),
                 std::inserter(groups, groups.begin()), [&root_grp](const auto group_name) {
                   const auto name = std::string{group_name};
                   auto group_obj = root_grp().createGroup(std::string{group_name});

                   return std::make_pair(name, Group{root_grp, group_obj});
                 });
  return groups;
}

template <typename PixelT>
inline auto File::create_datasets(RootGroup &root_grp, const Reference &chroms,
                                  std::size_t cache_size_bytes, double w0) -> DatasetMap {
  DatasetMap datasets(MANDATORY_DATASET_NAMES.size() + 1);

  const std::size_t num_pixel_datasets = 3;
  const std::size_t num_read_once_dataset = MANDATORY_DATASET_NAMES.size() - num_pixel_datasets;

  const std::size_t read_once_cache_size = DEFAULT_HDF5_DATASET_CACHE_SIZE;
  const std::size_t pixel_dataset_cache_size =
      (cache_size_bytes - (read_once_cache_size * num_read_once_dataset)) / num_pixel_datasets;

  const auto default_aprop =
      Dataset::init_access_props(DEFAULT_HDF5_CHUNK_SIZE, read_once_cache_size, 1.0);
  const auto pixels_aprop = Dataset::init_access_props(
      DEFAULT_HDF5_CHUNK_SIZE, (std::max(read_once_cache_size, pixel_dataset_cache_size)), w0);

  auto create_dataset = [&](const auto &path, const auto &type, auto aprop) {
    using T = remove_cvref_t<decltype(type)>;
    if constexpr (is_string_v<T>) {
      const auto &chrom_with_longest_name = chroms.chromosome_with_longest_name();
      datasets.emplace(path, Dataset{root_grp, path, chrom_with_longest_name.name(),
                                     HighFive::DataSpace::UNLIMITED, aprop});
    } else {
      datasets.emplace(path, Dataset{root_grp, path, type, HighFive::DataSpace::UNLIMITED, aprop});
    }
  };

  create_dataset("chroms/name", std::string{}, default_aprop);
  create_dataset("chroms/length", std::int32_t{}, default_aprop);

  create_dataset("bins/chrom", std::int32_t{}, default_aprop);
  create_dataset("bins/start", std::int32_t{}, default_aprop);
  create_dataset("bins/end", std::int32_t{}, default_aprop);

  create_dataset("pixels/bin1_id", std::int64_t{}, pixels_aprop);
  create_dataset("pixels/bin2_id", std::int64_t{}, pixels_aprop);
  create_dataset("pixels/count", PixelT{}, pixels_aprop);

  create_dataset("indexes/bin1_offset", std::int64_t{}, default_aprop);
  create_dataset("indexes/chrom_offset", std::int64_t{}, default_aprop);

  assert(datasets.size() == MANDATORY_DATASET_NAMES.size());

  return datasets;
}

inline void File::write_standard_attributes(RootGroup &root_grp,
                                            const StandardAttributes &attributes,
                                            bool skip_sentinel_attr) {
  assert(attributes.bin_size != 0);
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (attributes.assembly) {
    Attribute::write(root_grp(), "assembly", *attributes.assembly);
  }
  Attribute::write(root_grp(), "bin-size", attributes.bin_size);
  Attribute::write(root_grp(), "bin-type", *attributes.bin_type);            // NOLINT
  Attribute::write(root_grp(), "creation-date", *attributes.creation_date);  // NOLINT
  Attribute::write(root_grp(), "format", std::string{COOL_MAGIC});
  Attribute::write(root_grp(), "format-url", *attributes.format_url);  // NOLINT
  if (!skip_sentinel_attr) {
    static_assert(internal::SENTINEL_ATTR_NAME == "format-version");
    Attribute::write(root_grp(), "format-version", attributes.format_version);
  }
  Attribute::write(root_grp(), "generated-by", *attributes.generated_by);  // NOLINT
  Attribute::write(root_grp(), "metadata", *attributes.metadata);          // NOLINT
  Attribute::write(root_grp(), "nbins", *attributes.nbins);                // NOLINT
  Attribute::write(root_grp(), "nchroms", *attributes.nchroms);            // NOLINT
  Attribute::write(root_grp(), "nnz", *attributes.nnz);                    // NOLINT
  Attribute::write(root_grp(), "storage-mode", *attributes.storage_mode);  // NOLINT
  assert(attributes.sum.has_value());
  assert(attributes.cis.has_value());
  std::visit([&](const auto sum) { Attribute::write(root_grp(), "sum", sum); }, *attributes.sum);
  std::visit([&](const auto cis) { Attribute::write(root_grp(), "cis", cis); }, *attributes.cis);
}

inline void File::write_attributes(bool skip_sentinel_attr) {
  assert(this->_attrs.nbins == this->bins().size());
  assert(this->_attrs.nchroms == this->chromosomes().size());
  assert(this->_attrs.nnz == this->_datasets.at("pixels/count").size());

  File::write_standard_attributes(this->_root_group, this->_attrs, skip_sentinel_attr);
  this->flush();
  if (skip_sentinel_attr) {
    using T [[maybe_unused]] = remove_cvref_t<decltype(internal::SENTINEL_ATTR_VALUE)>;
    assert(Attribute::read<T>(this->_root_group(), internal::SENTINEL_ATTR_NAME) ==
           internal::SENTINEL_ATTR_VALUE);
    Attribute::write(this->_root_group(), "format-version", this->_attrs.format_version, true);
    this->flush();
  }
}

inline void File::write_chromosomes() {
  assert(this->_datasets.contains("chroms/name"));
  assert(this->_datasets.contains("chroms/length"));
  assert(!this->chromosomes().empty());

  File::write_chromosomes(this->dataset("chroms/name"), this->dataset("chroms/length"),
                          this->chromosomes().begin(), this->chromosomes().end());

  this->_attrs.nchroms = this->chromosomes().size();
}

template <typename ChromIt, typename UnaryOperation, typename>
inline void File::write_chromosomes(Dataset &name_dset, Dataset &size_dset, ChromIt first_chrom,
                                    ChromIt last_chrom, UnaryOperation op) {
  const auto num_chroms = std::distance(first_chrom, last_chrom);
  if (num_chroms == 0) {
    return;
  }

  try {
    name_dset.write(first_chrom, last_chrom, 0, true,
                    [&](const auto &chrom) { return std::string{op(chrom).name()}; });
  } catch (const HighFive::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write {} chromosome name(s) to \"{}\": {}"), num_chroms,
                    name_dset.uri(), e.what()));
  }
  try {
    size_dset.write(first_chrom, last_chrom, 0, true,
                    [&](const auto &chrom) { return op(chrom).size(); });
  } catch (const HighFive::Exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Failed to write {} chromosome size(s) to \"{}\": {}"), num_chroms,
                    size_dset.uri(), e.what()));
  }

  assert(name_dset.size() == static_cast<std::size_t>(num_chroms));
  assert(size_dset.size() == static_cast<std::size_t>(num_chroms));
}

inline void File::write_bin_table() {
  File::write_bin_table(this->dataset("bins/chrom"), this->dataset("bins/start"),
                        this->dataset("bins/end"), this->bins());

  this->_attrs.nbins = this->bins().size();
}

inline void File::write_bin_table(Dataset &chrom_dset, Dataset &start_dset, Dataset &end_dset,
                                  const BinTable &bin_table) {
  assert(!bin_table.empty());

  chrom_dset.write(bin_table.begin(), bin_table.end(), 0, true,
                   [&](const Bin &bin) { return bin.chrom().id(); });

  start_dset.write(bin_table.begin(), bin_table.end(), 0, true,
                   [&](const Bin &bin) { return bin.start(); });

  end_dset.write(bin_table.begin(), bin_table.end(), 0, true,
                 [&](const Bin &bin) { return bin.end(); });

  assert(chrom_dset.size() == bin_table.size());
  assert(start_dset.size() == bin_table.size());
  assert(end_dset.size() == bin_table.size());
}

template <typename PixelIt>
inline void File::update_indexes(PixelIt first_pixel, PixelIt last_pixel) {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = decltype(std::declval<PixelT>().count);

  if (first_pixel == last_pixel) {
    return;
  }

  auto nnz = static_cast<std::uint64_t>(*this->_attrs.nnz);
  PixelCoordinates first_pixel_in_row(this->get_last_bin_written());

  std::for_each(first_pixel, last_pixel, [&](const Pixel<T> &p) {
    if (first_pixel_in_row.bin1.start() != p.coords.bin1.start()) {
      first_pixel_in_row = p.coords;
      this->index().set_offset_by_bin_id(first_pixel_in_row.bin1.id(), nnz);
    }
    nnz++;
  });
}

inline void File::write_indexes() {
  assert(this->_attrs.nnz.has_value());
  this->index().finalize(static_cast<std::uint64_t>(*this->_attrs.nnz));
  File::write_indexes(this->dataset("indexes/chrom_offset"), this->dataset("indexes/bin1_offset"),
                      this->index());
}

inline void File::write_indexes(Dataset &chrom_offset_dset, Dataset &bin_offset_dset,
                                const Index &idx) {
  chrom_offset_dset.write(idx.compute_chrom_offsets(), 0, true);

  bin_offset_dset.write(idx.begin(), idx.end(), 0, true);

  assert(chrom_offset_dset.size() == idx.num_chromosomes() + 1);
  assert(bin_offset_dset.size() == idx.size() + 1);
}

inline void File::write_sentinel_attr(HighFive::Group grp) {
  assert(!check_sentinel_attr(grp));

  Attribute::write(grp, internal::SENTINEL_ATTR_NAME, internal::SENTINEL_ATTR_VALUE, true);
  grp.getFile().flush();
}

inline void File::write_sentinel_attr() { File::write_sentinel_attr(this->_root_group()); }

}  // namespace hictk
