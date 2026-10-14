// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/common.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/reference.hpp"

namespace hictk::cooler {

template <typename PixelT>
inline File::File(RootGroup entrypoint, BinTable bins, [[maybe_unused]] PixelT pixel,
                  Attributes attributes, std::size_t cache_size_bytes,
                  std::uint32_t compression_lvl, double w0)
    : _mode(HighFive::File::ReadWrite),
      _root_group(std::move(entrypoint)),
      _groups(create_groups(_root_group)),
      _datasets(create_datasets<PixelT>(_root_group, bins.chromosomes(), cache_size_bytes,
                                        compression_lvl, w0)),
      _attrs(std::move(attributes)),
      _pixel_variant(PixelT(0)),
      _bins(std::make_shared<const BinTable>(std::move(bins))),
      _index(std::make_shared<Index>(_bins)),
      _finalize(true) {
  assert(!_bins->empty());
  assert(!chromosomes().empty());
  assert(!_index->empty());
  assert(std::holds_alternative<PixelT>(_pixel_variant));
  write_chromosomes();
  write_bin_table();

  if constexpr (std::is_floating_point_v<PixelT>) {
    _attrs.sum = 0.0;
    _attrs.cis = 0.0;
  } else {
    _attrs.sum = std::int64_t{0};
    _attrs.cis = std::int64_t{0};
  }

  write_sentinel_attr();
}

template <typename PixelT>
inline File::File(RootGroup entrypoint, [[maybe_unused]] PixelT pixel, Attributes attributes,
                  std::size_t cache_size_bytes, double w0)
    : _mode(HighFive::File::ReadWrite),
      _root_group(std::move(entrypoint)),
      _attrs(std::move(attributes)),
      _pixel_variant(PixelT(0)),
      _finalize(true) {
  _groups = open_groups(_root_group);
  _datasets = open_datasets(_root_group, cache_size_bytes, w0);

  if (std::is_floating_point_v<PixelT>) {
    _attrs.sum = 0.0;
    _attrs.cis = 0.0;
  } else {
    _attrs.sum = std::int64_t{0};
    _attrs.cis = std::int64_t{0};
  }

  _bins = std::make_shared<BinTable>(init_bin_table(_datasets, _attrs.bin_type, _attrs.bin_size));
  _index = std::make_shared<Index>(_bins);

  assert(std::holds_alternative<PixelT>(_pixel_variant));
  assert(!_bins->empty());
  assert(!chromosomes().empty());
  assert(!_index->empty());

  write_sentinel_attr();
}

template <typename PixelT>
inline File File::create(std::string_view uri, const Reference &chroms, std::uint32_t bin_size,
                         bool overwrite_if_exists, Attributes attributes,
                         std::size_t cache_size_bytes, std::uint32_t compression_lvl) {
  return File::create<PixelT>(uri, BinTable(chroms, bin_size), overwrite_if_exists,
                              std::move(attributes), cache_size_bytes, compression_lvl);
}

template <typename PixelT>
inline File File::create(std::string_view uri, BinTable bins, bool overwrite_if_exists,
                         Attributes attributes, std::size_t cache_size_bytes,
                         std::uint32_t compression_lvl) {
  try {
    const auto [file_path, root_path] = parse_cooler_uri(uri);
    const auto uri_is_file_path = root_path.empty() || root_path == "/";

    // URI is like myfile.mcool::/resolutions/100, but myfile.mcool does not exist
    if (!uri_is_file_path && !std::filesystem::exists(file_path)) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("parent file \"{}\" does not exist.\n"
                     "Did you forget to create the parent file with e.g. init_mcool()?"),
          uri, file_path));
    }

    // URI points to an existing file, but overwrite_if_exists=false
    if (!overwrite_if_exists && uri_is_file_path && std::filesystem::exists(file_path)) {
      throw std::runtime_error("URI points to an existing file");
    }

    auto mode = overwrite_if_exists ? HighFive::File::Overwrite : HighFive::File::Create;

    // File exists but cooler may not
    if (std::filesystem::exists(file_path) && !uri_is_file_path) {
      mode = HighFive::File::ReadWrite;
    }

    {
      auto fp = open_file(uri, mode, false);
      auto root_group = open_or_create_root_group(fp, uri);
      if (!uri_is_file_path && utils::is_cooler(root_group())) {
        if (overwrite_if_exists) {
          throw std::runtime_error(fmt::format(
              FMT_STRING("overwriting cooler nested inside .mcool or .scool is not yet supported.\n"
                         "Path to parent file: \"{}\"\""
                         "Path to nested cooler: \"{}\""),
              file_path, root_path));
        }
      }
      assert(!utils::is_cooler(root_group()));
    }

    return create<PixelT>(
        open_or_create_root_group(open_file(uri, HighFive::File::ReadWrite, false), uri),
        std::move(bins), std::move(attributes), cache_size_bytes, compression_lvl);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Cannot create cooler at the following URI: \"{}\". Reason: {}"),
                    uri, e.what()));
  }
}

template <typename PixelT>
inline File File::create(RootGroup entrypoint, const Reference &chroms, std::uint32_t bin_size,
                         Attributes attributes, std::size_t cache_size_bytes,
                         std::uint32_t compression_lvl) {
  return File::create<PixelT>(std::move(entrypoint), BinTable(chroms, bin_size),
                              std::move(attributes), cache_size_bytes, compression_lvl);
}

template <typename PixelT>
inline File File::create(RootGroup entrypoint, BinTable bins, Attributes attributes,
                         std::size_t cache_size_bytes, std::uint32_t compression_lvl) {
  static_assert(std::is_arithmetic_v<PixelT>);
  attributes.bin_type = bins.type();
  attributes.bin_size = bins.resolution();

  try {
    if (utils::is_cooler(entrypoint())) {
      throw std::runtime_error("URI points to an already existing cooler.");
    }
    return File(entrypoint, std::move(bins), PixelT(0), attributes, cache_size_bytes,
                compression_lvl, DEFAULT_HDF5_CACHE_W0);

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Cannot create cooler at the following URI: \"{}\". Reason: {}"),
                    entrypoint.uri(), e.what()));
  }
}

constexpr File::operator bool() const noexcept { return bool(_root_group); }

constexpr bool File::operator!() const noexcept { return !bool(*this); }

template <typename N, bool cis>
inline void File::update_pixel_sum(N partial_sum) {
  static_assert(std::is_arithmetic_v<N>);

  auto &buff = cis ? _attrs.cis : _attrs.sum;
  assert(buff.has_value());
  if constexpr (std::is_floating_point_v<N>) {
    std::get<double>(*buff) += conditional_static_cast<double>(partial_sum);
  } else {
    assert(std::is_integral_v<N>);
    std::get<std::int64_t>(*buff) += conditional_static_cast<std::int64_t>(partial_sum);
  }
}

}  // namespace hictk::cooler
