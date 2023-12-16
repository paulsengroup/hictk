// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

template <typename PixelIt>
void File::append_bins(Dataset &bin1_dset, Dataset &bin2_dset, PixelIt first_pixel,
                       PixelIt last_pixel) {
  using PixelT = std::decay_t<decltype(*first_pixel)>;
  using T = remove_cvref_t<decltype(first_pixel->count)>;

  if constexpr (std::is_same_v<PixelT, Pixel<T>>) {
    bin1_dset.append(first_pixel, last_pixel,
                     [&](const auto &pixel) { return pixel.coords.bin1.id(); });

    bin2_dset.append(first_pixel, last_pixel,
                     [&](const auto &pixel) { return pixel.coords.bin2.id(); });
  } else {
    bin1_dset.append(first_pixel, last_pixel, [&](const auto &pixel) { return pixel.bin1_id; });

    bin2_dset.append(first_pixel, last_pixel, [&](const auto &pixel) { return pixel.bin2_id; });
  }
}

namespace internal {
template <typename PixelT, typename T>
PixelCoordinates pixel_to_coords(const PixelT &pixel, const BinTable &bins) {
  if constexpr (std::is_same_v<PixelT, Pixel<T>>) {
    return pixel.coords;
  } else {
    return PixelCoordinates{bins.at(pixel.bin1_id), bins.at(pixel.bin2_id)};
  }
}
}  // namespace internal

template <typename PixelIt, typename N>
void File::append_counts(Dataset &dset, const BinTable &bins, PixelIt first_pixel,
                         PixelIt last_pixel, N &sum, N &cis_sum) {
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = remove_cvref_t<decltype(first_pixel->count)>;

  sum = 0;
  cis_sum = 0;

  dset.append(first_pixel, last_pixel, [&](const auto &pixel) {
    if (pixel.count == 0) {
      // Defining pixel_to_coords as a lambda and calling fmt::format to format the error message
      // causes a compiler segfault when using conda-forge's compiler toolchain
      throw std::runtime_error("Found pixel with 0 interactions: " +
                               fmt::to_string(internal::pixel_to_coords<PixelT, T>(pixel, bins)));
    }

    sum += pixel.count;
    if constexpr (std::is_same_v<PixelT, Pixel<T>>) {
      if (pixel.coords.bin1.chrom().id() == pixel.coords.bin2.chrom().id()) {
        cis_sum += pixel.count;
      }
    } else {
      const auto chrom1 = bins.at(pixel.bin1_id).chrom();
      const auto chrom2 = bins.at(pixel.bin2_id).chrom();
      if (chrom1 == chrom2) {
        cis_sum += pixel.count;
      }
    }
    return pixel.count;
  });
}

template <typename PixelIt, typename>
inline void File::append_pixels(PixelIt first_pixel, PixelIt last_pixel, bool validate) {
  using PixelT = std::decay_t<decltype(*first_pixel)>;
  using T = remove_cvref_t<decltype(first_pixel->count)>;
  if constexpr (ndebug_not_defined()) {
    validate_pixel_type<T>();
  }

  update_indexes(first_pixel, last_pixel);

  if (validate) {
    if constexpr (std::is_same_v<PixelT, Pixel<T>>) {
      validate_pixels_before_append(first_pixel, last_pixel);
    } else {
      validate_thin_pixels_before_append(first_pixel, last_pixel);
    }
  }

  File::append_bins(dataset("pixels/bin1_id"), dataset("pixels/bin2_id"), first_pixel, last_pixel);

  // NOLINTBEGIN(*-avoid-non-const-global-variables)
  T sum{};
  T cis_sum{};
  // NOLINTEND(*-avoid-non-const-global-variables)
  File::append_counts(dataset("pixels/count"), bins(), first_pixel, last_pixel, sum, cis_sum);
  _attrs.nnz = dataset("pixels/bin1_id").size();

  update_pixel_sum(sum);
  update_pixel_sum<T, true>(cis_sum);
}

inline void File::flush() { _root_group().getFile().flush(); }

template <typename It>
inline void File::write_weights(std::string_view uri, std::string_view name, It first_weight,
                                It last_weight, bool overwrite_if_exists, bool divisive) {
  File(open_or_create_root_group(open_file(uri, HighFive::File::ReadWrite, true), uri),
       HighFive::File::ReadWrite, DEFAULT_HDF5_CACHE_SIZE, DEFAULT_HDF5_CACHE_W0, true)
      .write_weights(name, first_weight, last_weight, overwrite_if_exists, divisive);
}

template <typename It>
inline void File::write_weights(std::string_view name, It first_weight, It last_weight,
                                bool overwrite_if_exists, bool divisive) {
  if (name.empty()) {
    throw std::runtime_error("weight name is empty");
  }

  if (_mode == HighFive::File::ReadOnly) {
    throw std::runtime_error("File::write_weights() was called on a file open in read-only mode");
  }

  const auto weights_shape = static_cast<std::size_t>(std::distance(first_weight, last_weight));
  const auto expected_weights_shape = bins().size();
  if (weights_shape != expected_weights_shape) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Invalid weight shape, expected {} values, found {}"),
                    expected_weights_shape, weights_shape));
  }

  const auto path = fmt::format(FMT_STRING("bins/{}"), name);
  const auto existing = _root_group().exist(path);
  if (!overwrite_if_exists && existing) {
    throw std::runtime_error(fmt::format(FMT_STRING("dataset \"{}\" already exists"), path));
  }

  typename std::iterator_traits<It>::value_type buff{};
  auto dset = existing ? Dataset(_root_group, _root_group().getDataSet(path))
                       : Dataset(_root_group, path, buff, HighFive::DataSpace::UNLIMITED);

  dset.resize(weights_shape);
  if (weights_shape != 0) {
    dset.write(first_weight, last_weight);
  }

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

inline auto File::create_groups(RootGroup &root_grp, Group chroms_grp, Group bins_grp) -> GroupMap {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  GroupMap groups(MANDATORY_GROUP_NAMES.size() + 1);
  root_grp().createHardLink("chroms", chroms_grp());
  root_grp().createHardLink("bins", bins_grp());

  groups.emplace(root_grp.hdf5_path(), Group{root_grp, root_grp()});
  groups.emplace("chroms", Group{root_grp, root_grp().getGroup("chroms")});
  groups.emplace("bins", Group{root_grp, root_grp().getGroup("bins")});

  groups.emplace("pixels", Group{root_grp, root_grp().createGroup("pixels")});
  groups.emplace("indexes", Group{root_grp, root_grp().createGroup("indexes")});

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
      DEFAULT_HDF5_CHUNK_SIZE, ((std::max)(read_once_cache_size, pixel_dataset_cache_size)), w0);

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

inline void File::write_standard_attributes(RootGroup &root_grp, const Attributes &attributes,
                                            bool skip_sentinel_attr) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (attributes.assembly) {
    Attribute::write(root_grp(), "assembly", *attributes.assembly);
  }
  if (attributes.bin_size == 0) {
    assert(attributes.bin_type.value() == "variable");
    Attribute::write(root_grp(), "bin-size", "null");
  } else {
    assert(attributes.bin_type.value() == "fixed");
    Attribute::write(root_grp(), "bin-size", attributes.bin_size);
  }
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
  assert(_attrs.nbins == bins().size());
  assert(_attrs.nchroms == chromosomes().size());
  assert(_attrs.nnz == _datasets.at("pixels/count").size());

  File::write_standard_attributes(_root_group, _attrs, skip_sentinel_attr);
  flush();
  if (skip_sentinel_attr) {
    using T [[maybe_unused]] = remove_cvref_t<decltype(internal::SENTINEL_ATTR_VALUE)>;
    assert(Attribute::read<T>(_root_group(), internal::SENTINEL_ATTR_NAME) ==
           internal::SENTINEL_ATTR_VALUE);
    Attribute::write(_root_group(), "format-version", _attrs.format_version, true);
    flush();
  }
}

inline void File::write_chromosomes() {
  assert(_datasets.contains("chroms/name"));
  assert(_datasets.contains("chroms/length"));
  assert(!chromosomes().empty());

  File::write_chromosomes(dataset("chroms/name"), dataset("chroms/length"), chromosomes().begin(),
                          chromosomes().end());

  _attrs.nchroms = static_cast<std::int32_t>(chromosomes().size());
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
  File::write_bin_table(dataset("bins/chrom"), dataset("bins/start"), dataset("bins/end"), bins());

  _attrs.nbins = bins().size();
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
  using PixelT = std::decay_t<decltype(*first_pixel)>;
  using T = remove_cvref_t<decltype(first_pixel->count)>;

  if (first_pixel == last_pixel) {
    return;
  }

  auto nnz = static_cast<std::uint64_t>(*_attrs.nnz);
  PixelCoordinates first_pixel_in_row(get_last_bin_written());

  if constexpr (std::is_same_v<PixelT, Pixel<T>>) {
    std::for_each(first_pixel, last_pixel, [&](const Pixel<T> &p) {
      if (first_pixel_in_row.bin1 != p.coords.bin1) {
        first_pixel_in_row = p.coords;
        index().set_offset_by_bin_id(first_pixel_in_row.bin1.id(), nnz);
      }
      nnz++;
    });
  } else {
    std::for_each(first_pixel, last_pixel, [&](const ThinPixel<T> &p) {
      if (first_pixel_in_row.bin1.id() != p.bin1_id) {
        first_pixel_in_row = {bins().at(p.bin1_id), bins().at(p.bin2_id)};
        index().set_offset_by_bin_id(first_pixel_in_row.bin1.id(), nnz);
      }
      nnz++;
    });
  }
}

inline void File::write_indexes() {
  assert(_attrs.nnz.has_value());
  index().finalize(static_cast<std::uint64_t>(*_attrs.nnz));
  File::write_indexes(dataset("indexes/chrom_offset"), dataset("indexes/bin1_offset"), index());
}

inline void File::write_indexes(Dataset &chrom_offset_dset, Dataset &bin_offset_dset,
                                const Index &idx) {
  chrom_offset_dset.write(idx.compute_chrom_offsets(), 0, true);

  bin_offset_dset.write(idx.begin(), idx.end(), 0, true);

  assert(chrom_offset_dset.size() == idx.chromosomes().size() + 1);
  assert(bin_offset_dset.size() == idx.size() + 1);
}

inline void File::write_sentinel_attr(HighFive::Group grp) {
  assert(!check_sentinel_attr(grp));

  Attribute::write(grp, internal::SENTINEL_ATTR_NAME, internal::SENTINEL_ATTR_VALUE, true);
  grp.getFile().flush();
}

inline void File::write_sentinel_attr() { File::write_sentinel_attr(_root_group()); }

}  // namespace hictk::cooler
