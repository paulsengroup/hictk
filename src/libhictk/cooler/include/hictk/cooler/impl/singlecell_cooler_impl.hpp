// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <highfive/H5DataSpace.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/cooler/utils.hpp"

namespace hictk::cooler {

constexpr const phmap::btree_set<std::string>& SingleCellFile::cells() const noexcept {
  return _cells;
}
constexpr const SingleCellAttributes& SingleCellFile::attributes() const noexcept { return _attrs; }

template <typename N>
inline File SingleCellFile::create_cell(std::string_view cell, Attributes attrs,
                                        std::size_t cache_size_bytes,
                                        std::uint32_t compression_lvl) {
  assert(_root_grp);
  if (_cells.contains(cell)) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("failed to create cell \"{}\" in file {}: cell already exists"), cell, path()));
  }

  auto sclr_attrs = read_standard_attributes((*_root_grp)().getFile(), true);

  attrs.bin_size = sclr_attrs.bin_size;
  attrs.bin_type = sclr_attrs.bin_type;
  attrs.format = sclr_attrs.format;
  attrs.format_version = sclr_attrs.format_version;
  attrs.creation_date = sclr_attrs.creation_date;
  attrs.generated_by = sclr_attrs.generated_by;
  attrs.assembly = sclr_attrs.assembly;
  attrs.metadata = sclr_attrs.metadata;
  attrs.format_url = sclr_attrs.format_url;
  attrs.nbins = sclr_attrs.nbins;
  attrs.nchroms = sclr_attrs.nchroms;
  attrs.storage_mode = sclr_attrs.storage_mode;

  _cells.emplace(cell);
  Attribute::write((*_root_grp)(), "ncells", static_cast<std::int64_t>(_cells.size()), true);

  RootGroup entrypoint{
      (*_root_grp)().createGroup(fmt::format(FMT_STRING("cells/{}"), std::string{cell}))};

  std::ignore =
      File::create_groups(entrypoint, Group{*_root_grp, (*_root_grp)().getGroup("/chroms")},
                          Group{*_root_grp, (*_root_grp)().getGroup("/bins")});

  create_cell_datasets<N>(entrypoint, cache_size_bytes, compression_lvl, DEFAULT_HDF5_CACHE_W0);

  return {entrypoint, N{}, std::move(attrs), cache_size_bytes, DEFAULT_HDF5_CACHE_W0};
}

template <typename N>
inline File SingleCellFile::aggregate(std::string_view uri, bool overwrite_if_exists,
                                      std::uint32_t compression_lvl, std::size_t chunk_size,
                                      std::size_t update_frequency) const {
  if (_cells.size() == 1) {
    if (overwrite_if_exists && std::filesystem::exists(uri)) {
      std::filesystem::remove(uri);
    }
    utils::copy(open(*_cells.begin()).uri(), uri);
    return File(uri);
  }

  std::vector<cooler::PixelSelector::iterator<N>> heads{};
  std::vector<cooler::PixelSelector::iterator<N>> tails{};

  std::for_each(_cells.begin(), _cells.end(), [&](const auto& cell) {
    auto clr = open(cell);
    auto first = clr.template begin<N>();
    auto last = clr.template end<N>();
    if (first != last) {
      heads.emplace_back(std::move(first));
      tails.emplace_back(std::move(last));
    }
  });

  utils::merge(heads, tails, bins(), uri, attributes().assembly.value_or("unknown"),
               overwrite_if_exists, chunk_size, update_frequency, compression_lvl);

  return File(uri);
}

template <typename PixelT>
inline void SingleCellFile::create_cell_datasets(RootGroup& root_grp, std::size_t cache_size_bytes,
                                                 std::uint32_t compression_lvl, double w0) {
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

  auto create_dataset = [&](const auto& path, const auto& type, const auto& aprop,
                            const auto& cprop) {
    Dataset{root_grp, path, type, HighFive::DataSpace::UNLIMITED, aprop, cprop};
  };

  create_dataset("pixels/bin1_id", std::int64_t{}, pixels_aprop, default_cprop);
  create_dataset("pixels/bin2_id", std::int64_t{}, pixels_aprop, default_cprop);
  create_dataset("pixels/count", PixelT{}, pixels_aprop, default_cprop);

  create_dataset("indexes/bin1_offset", std::int64_t{}, default_aprop, default_cprop);
  create_dataset("indexes/chrom_offset", std::int64_t{}, default_aprop, default_cprop);
}

}  // namespace hictk::cooler
