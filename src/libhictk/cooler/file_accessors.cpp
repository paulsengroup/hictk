// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "hictk/cooler/cooler.hpp"
#include "hictk/suppress_warnings.hpp"
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
HICTK_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/numeric_variant.hpp"
#include "hictk/reference.hpp"

namespace hictk::cooler {

std::string File::uri() const {
  if (!*this) {
    return "";
  }
  if (hdf5_path() == "/") {
    return path();
  }
  return fmt::format(FMT_STRING("{}::{}"), path(), hdf5_path());
}

std::string File::hdf5_path() const { return _root_group.hdf5_path(); }
std::string File::path() const {
  if (!*this) {
    return "";
  }
  return _root_group().getFile().getName();
}

auto File::chromosomes() const noexcept -> const Reference & { return bins().chromosomes(); }

auto File::bins() const noexcept -> const BinTable & {
  assert(_bins);
  return *_bins;
}

auto File::bins_ptr() const noexcept -> std::shared_ptr<const BinTable> { return _bins; }

std::uint32_t File::resolution() const noexcept { return _attrs.bin_size; }
std::uint64_t File::nbins() const { return bins().size(); }
std::uint64_t File::nchroms() const { return chromosomes().size(); }
std::uint64_t File::nnz() const { return dataset("pixels/count").size(); }

auto File::attributes() const noexcept -> const Attributes & { return _attrs; }

HighFive::File File::file_handle() { return _root_group().getFile(); }

const HighFive::File &File::file_handle() const { return _root_group().getFile(); }

auto File::group(std::string_view group_name) -> Group & {
  try {
    return _groups.at(std::string{group_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Group \"{}\" does not exists!"), group_name));
  }
}

auto File::group(std::string_view group_name) const -> const Group & {
  try {
    return _groups.at(std::string{group_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("Group \"{}\" does not exists!"), group_name));
  }
}

auto File::dataset(std::string_view dataset_name) -> Dataset & {
  try {
    if (dataset_name.front() == '/') {
      dataset_name = dataset_name.substr(1);
    }
    return _datasets.at(std::string{dataset_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Dataset \"{}\" does not exists!"), dataset_name));
  }
}

auto File::dataset(std::string_view dataset_name) const -> const Dataset & {
  try {
    if (dataset_name.front() == '/') {
      dataset_name = dataset_name.substr(1);
    }
    return _datasets.at(std::string{dataset_name});
  } catch ([[maybe_unused]] const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Dataset \"{}\" does not exists!"), dataset_name));
  }
}

bool File::has_normalization(const balancing::Method &normalization) const {
  const auto dset_path =
      fmt::format(FMT_STRING("{}/{}"), _groups.at("bins")().getPath(), normalization.to_string());
  if (_weights.contains(dset_path)) {
    return true;
  }

  return _root_group().exist(dset_path);
}

std::vector<balancing::Method> File::avail_normalizations() const {
  const phmap::flat_hash_set<std::string> bin_table_dsets{"chrom", "start", "end"};

  std::vector<balancing::Method> norms{};
  for (const auto &dset : group("bins")().listObjectNames(HighFive::IndexType::NAME)) {
    if (bin_table_dsets.contains(dset)) {
      continue;
    }

    norms.emplace_back(dset);
  }

  return norms;
}

const hictk::internal::NumericVariant &File::pixel_variant() const noexcept {
  return _pixel_variant;
}

bool File::has_signed_pixels() const noexcept {
  // clang-format off
  return has_pixel_of_type<std::int8_t>()  ||
         has_pixel_of_type<std::int16_t>() ||
         has_pixel_of_type<std::int32_t>() ||
         has_pixel_of_type<std::int64_t>();
  // clang-format on
}

bool File::has_unsigned_pixels() const noexcept {
  // clang-format off
  return has_pixel_of_type<std::uint8_t>()  ||
         has_pixel_of_type<std::uint16_t>() ||
         has_pixel_of_type<std::uint32_t>() ||
         has_pixel_of_type<std::uint64_t>();
  // clang-format on
}

bool File::has_integral_pixels() const noexcept {
  return has_signed_pixels() || has_unsigned_pixels();
}

bool File::has_float_pixels() const noexcept {
  // clang-format off
  return has_pixel_of_type<float>()  ||
         has_pixel_of_type<double>() ||
         has_pixel_of_type<long double>();
  // clang-format on
}

auto File::index() noexcept -> Index & {
  assert(_index);
  return *_index;
}

auto File::index() const noexcept -> const Index & {
  assert(_index);
  return *_index;
}

}  // namespace hictk::cooler
