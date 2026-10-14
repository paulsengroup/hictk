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
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

template <typename PixelIt>
void File::append_bins(Dataset &bin1_dset, Dataset &bin2_dset, const PixelIt &first_pixel,
                       const PixelIt &last_pixel) {
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

// NOLINTBEGIN(*-unnecessary-value-param)
template <typename PixelIt, typename N>
void File::append_counts(Dataset &dset, const BinTable &bins, PixelIt first_pixel,
                         PixelIt last_pixel, N &sum, N &cis_sum) {
  // NOLINTEND(*-unnecessary-value-param)
  using PixelT = typename std::iterator_traits<PixelIt>::value_type;
  using T = remove_cvref_t<decltype(first_pixel->count)>;

  sum = 0;
  cis_sum = 0;

  dset.append(std::move(first_pixel), std::move(last_pixel), [&](const auto &pixel) {
    if (pixel.count == 0) {
      // Defining pixel_to_coords as a lambda and calling fmt::format to format the error message
      // causes a compiler segfault when using conda-forge's compiler toolchain
      throw std::runtime_error("Found pixel with 0 interactions: " +
                               fmt::to_string(internal::pixel_to_coords<PixelT, T>(pixel, bins)));
    }

    sum += static_cast<N>(pixel.count);
    if constexpr (std::is_same_v<PixelT, Pixel<T>>) {
      if (pixel.coords.bin1.chrom().id() == pixel.coords.bin2.chrom().id()) {
        cis_sum += static_cast<N>(pixel.count);
      }
    } else {
      const auto chrom1 = bins.at(pixel.bin1_id).chrom();
      const auto chrom2 = bins.at(pixel.bin2_id).chrom();
      if (chrom1 == chrom2) {
        cis_sum += static_cast<N>(pixel.count);
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

  if (validate) {
    if constexpr (std::is_same_v<PixelT, Pixel<T>>) {
      validate_pixels_before_append(first_pixel, last_pixel);
    } else {
      validate_thin_pixels_before_append(first_pixel, last_pixel);
    }
  }

  update_indexes(first_pixel, last_pixel);

  File::append_bins(dataset("pixels/bin1_id"), dataset("pixels/bin2_id"), first_pixel, last_pixel);

  // NOLINTBEGIN(*-avoid-non-const-global-variables)
  Attributes::SumVar sumv{};
  Attributes::SumVar cis_sumv{};
  if constexpr (std::is_floating_point_v<T>) {
    sumv = 0.0;
    cis_sumv = 0.0;
  } else {
    sumv = std::int64_t{0};
    cis_sumv = std::int64_t{0};
  }

  // NOLINTEND(*-avoid-non-const-global-variables)
  std::visit(
      [&](auto &sum) {
        using N = remove_cvref_t<decltype(sum)>;
        auto &cis_sum = std::get<N>(cis_sumv);

        File::append_counts(dataset("pixels/count"), bins(), std::move(first_pixel),
                            std::move(last_pixel), sum, cis_sum);
        _attrs.nnz = dataset("pixels/bin1_id").size();

        update_pixel_sum(sum);
        update_pixel_sum<N, true>(cis_sum);
      },
      sumv);
}

template <typename It>
inline void File::write_weights(std::string_view uri, std::string_view name, It first_weight,
                                It last_weight, bool overwrite_if_exists, bool divisive) {
  File(open_or_create_root_group(open_file(uri, HighFive::File::ReadWrite, true), uri),
       HighFive::File::ReadWrite, DEFAULT_HDF5_CACHE_SIZE * 4, DEFAULT_HDF5_CACHE_W0, true)
      .write_weights(name, first_weight, last_weight, overwrite_if_exists, divisive);
}

template <typename It>
inline void File::write_weights(std::string_view name, It first_weight, It last_weight,
                                bool overwrite_if_exists, bool divisive) {
  if (name.empty()) {
    throw std::runtime_error("weight name is empty");
  }

  if (name == "NONE") {
    throw std::runtime_error("caught attempt to write NONE weights");
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

  const typename std::iterator_traits<It>::value_type buff{};
  auto dset = existing ? Dataset(_root_group, _root_group().getDataSet(path))
                       : Dataset(_root_group, path, buff, HighFive::DataSpace::UNLIMITED);

  dset.resize(weights_shape);
  if (weights_shape != 0) {
    dset.write(first_weight, last_weight);
  }

  dset.write_attribute("divisive_weights", std::uint8_t{divisive}, overwrite_if_exists);
}

template <typename PixelT>
inline auto File::create_datasets(RootGroup &root_grp, const Reference &chroms,
                                  std::size_t cache_size_bytes, std::uint32_t compression_lvl,
                                  double w0) -> DatasetMap {
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

  const auto default_cprop = Dataset::init_create_props(compression_lvl, DEFAULT_HDF5_CHUNK_SIZE);

  auto create_dataset = [&](const auto &path, const auto &type, const auto &aprop,
                            const auto &cprop) {
    using T = remove_cvref_t<decltype(type)>;
    if constexpr (is_string_v<T>) {
      const auto &chrom_with_longest_name = chroms.chromosome_with_longest_name();
      datasets.emplace(path, Dataset{root_grp, path, chrom_with_longest_name.name(),
                                     HighFive::DataSpace::UNLIMITED, aprop, cprop});
    } else {
      datasets.emplace(path,
                       Dataset{root_grp, path, type, HighFive::DataSpace::UNLIMITED, aprop, cprop});
    }
  };

  create_dataset("chroms/name", std::string{}, default_aprop, default_cprop);
  create_dataset("chroms/length", std::int32_t{}, default_aprop, default_cprop);

  create_dataset("bins/chrom", std::int32_t{}, default_aprop, default_cprop);
  create_dataset("bins/start", std::int32_t{}, default_aprop, default_cprop);
  create_dataset("bins/end", std::int32_t{}, default_aprop, default_cprop);

  create_dataset("pixels/bin1_id", std::int64_t{}, pixels_aprop, default_cprop);
  create_dataset("pixels/bin2_id", std::int64_t{}, pixels_aprop, default_cprop);
  create_dataset("pixels/count", PixelT{}, pixels_aprop, default_cprop);

  create_dataset("indexes/bin1_offset", std::int64_t{}, default_aprop, default_cprop);
  create_dataset("indexes/chrom_offset", std::int64_t{}, default_aprop, default_cprop);

  assert(datasets.size() == MANDATORY_DATASET_NAMES.size());

  return datasets;
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

template <typename PixelIt>
inline void File::update_indexes(PixelIt first_pixel, PixelIt last_pixel) {
  using PixelT = std::decay_t<decltype(*first_pixel)>;
  using T = remove_cvref_t<decltype(first_pixel->count)>;

  if (first_pixel == last_pixel) {
    return;
  }

  assert(_attrs.nnz.has_value());  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  auto nnz = static_cast<std::uint64_t>(*_attrs.nnz);
  PixelCoordinates first_pixel_in_row(get_last_bin_written());

  if constexpr (std::is_same_v<PixelT, Pixel<T>>) {
    std::for_each(std::move(first_pixel), std::move(last_pixel), [&](const Pixel<T> &p) {
      if (first_pixel_in_row.bin1 != p.coords.bin1) {
        first_pixel_in_row = p.coords;
        index().set_offset_by_bin_id(first_pixel_in_row.bin1.id(), nnz);
      }
      nnz++;
    });
  } else {
    std::for_each(std::move(first_pixel), std::move(last_pixel), [&](const ThinPixel<T> &p) {
      if (first_pixel_in_row.bin1.id() != p.bin1_id) {
        first_pixel_in_row = {bins().at(p.bin1_id), bins().at(p.bin2_id)};
        index().set_offset_by_bin_id(first_pixel_in_row.bin1.id(), nnz);
      }
      nnz++;
    });
  }
}

}  // namespace hictk::cooler
