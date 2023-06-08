// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>
#include <memory>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/type_pretty_printer.hpp"
#include "hictk/variant_buff.hpp"

namespace hictk {

template <typename InputIt>
inline void init_mcool(std::string_view file_path, InputIt first_resolution,
                       InputIt last_resolution, bool force_overwrite) {
  using I = remove_cvref_t<decltype(*first_resolution)>;
  static_assert(std::is_integral_v<I>,
                "InputIt should be an iterator over a collection of integral numbers.");
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  const auto mode = force_overwrite ? HighFive::File::Truncate : HighFive::File::Create;
  HighFive::File fp(std::string{file_path}, mode);
  Attribute::write(fp, "format", std::string{MCOOL_MAGIC});
  Attribute::write(fp, "format-version", std::int64_t(3));

  auto res_group = fp.createGroup("/resolutions");
  std::for_each(first_resolution, last_resolution, [&](auto res) {
    assert(res > 0);
    auto cooler_root = res_group.createGroup(fmt::to_string(res));
  });
}

inline void init_mcool(std::string_view file_path, bool force_overwrite) {
  static constexpr std::array<std::uint64_t, 0> buff{};
  init_mcool(file_path, buff.begin(), buff.end(), force_overwrite);
}

// template <typename ChromSizeInputIt, typename CellIDInputIt>
// inline void init_scool(std::string_view file_path, ChromSizeInputIt first_chrom,
//                        ChromSizeInputIt last_chrom, CellIDInputIt first_cell_id,
//                        CellIDInputIt last_cell_id, std::uint32_t bin_size, bool force_overwrite)
//                        {
//   [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
//   const auto mode = force_overwrite ? IO_MODE::Truncate : IO_MODE::Create;
//   HighFive::File fp(std::string{file_path}, static_cast<unsigned>(mode));
//   fp.createAttribute("format", std::string{SCOOL_MAGIC});
//   fp.createAttribute("format-version", std::int64_t(3));
//
//   const auto bin_table = binnify(first_chrom, last_chrom, bin_size);
//
//   auto res_group = fp.createGroup("/resolutions");
//
// }

inline File::File(std::string_view uri, unsigned mode, std::size_t cache_size_bytes, double w0,
                  bool validate)
    : _mode(mode),
      _fp(std::make_unique<HighFive::File>(open_file(uri, _mode, validate))),
      _root_group(open_root_group(*_fp, uri)),
      _groups(open_groups(_root_group)),
      _datasets(open_datasets(_root_group, cache_size_bytes, w0)),
      _attrs(read_standard_attributes(_root_group)),
      _pixel_variant(detect_pixel_type(_root_group)),
      _bins(std::make_shared<BinTable>(
          import_chroms(_datasets.at("chroms/name"), _datasets.at("chroms/length"), false),
          this->bin_size())),
      _index(std::make_shared<Index>(
          import_indexes(_datasets.at("indexes/chrom_offset"), _datasets.at("indexes/bin1_offset"),
                         // NOLINTNEXTLINE
                         chromosomes(), _bins, static_cast<std::uint64_t>(*_attrs.nnz), false))) {
  assert(mode == HighFive::File::ReadOnly || mode == HighFive::File::ReadWrite);
  if (validate) {
    this->validate_bins();
  }
}

template <typename PixelT>
inline File::File(std::string_view uri, Reference chroms, [[maybe_unused]] PixelT pixel,
                  StandardAttributes attributes, std::size_t cache_size_bytes, double w0)
    : _mode(HighFive::File::ReadWrite),
      _fp(std::make_unique<HighFive::File>(open_file(uri, _mode, false))),
      _root_group(open_or_create_root_group(*_fp, uri)),
      _groups(create_groups(_root_group)),
      _datasets(create_datasets<PixelT>(_root_group, chroms, cache_size_bytes, w0)),
      _attrs(std::move(attributes)),
      _pixel_variant(PixelT(0)),
      _bins(std::make_shared<const BinTable>(std::move(chroms), this->bin_size())),
      _index(std::make_shared<Index>(_bins)),
      _finalize(true) {
  assert(this->bin_size() != 0);
  assert(!_bins->empty());
  assert(!chromosomes().empty());
  assert(!_index->empty());
  assert(std::holds_alternative<PixelT>(this->_pixel_variant));

  this->write_sentinel_attr();
}

inline File File::open_read_only(std::string_view uri, std::size_t cache_size_bytes,
                                 bool validate) {
  return File::open_read_only_random_access(uri, cache_size_bytes, validate);
}

inline File File::open_read_only_random_access(std::string_view uri, std::size_t cache_size_bytes,
                                               bool validate) {
  return File(uri, HighFive::File::ReadOnly, cache_size_bytes, DEFAULT_HDF5_CACHE_W0, validate);
}

inline File File::open_read_only_read_once(std::string_view uri, std::size_t cache_size_bytes,
                                           bool validate) {
  return File(uri, HighFive::File::ReadOnly, cache_size_bytes, 1.0, validate);
}

template <typename PixelT>
inline File File::create_new_cooler(std::string_view uri, const Reference &chroms,
                                    std::uint32_t bin_size, bool overwrite_if_exists,
                                    StandardAttributes attributes, std::size_t cache_size_bytes) {
  static_assert(std::is_arithmetic_v<PixelT>);
  if (bin_size == 0) {
    throw std::logic_error("bin_size cannot be zero.");
  }
  attributes.bin_size = bin_size;
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
    // At this point the parent file is guaranteed to exist, so we can always open it in ReadWrite
    // mode
    return File(uri, chroms, PixelT(0), attributes, cache_size_bytes);

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Cannot create cooler at the following URI: \"{}\". Reason: {}"),
                    uri, e.what()));
  }
}

inline File::~File() noexcept {
  try {
    this->finalize();
  } catch (const std::exception &e) {
    fmt::print(stderr, FMT_STRING("{}\n"), e.what());
  } catch (...) {
    fmt::print(stderr,
               FMT_STRING("An unknown error occurred while closing file {}. File is likely "
                          "corrupted or incomplete."),
               this->path());
  }
}

inline File::operator bool() const noexcept { return !!this->_fp; }

inline void File::open(std::string_view uri, bool validate) {
  *this = File::open_read_only(uri, validate);
}

template <typename PixelT>
inline void File::create(std::string_view uri, const hictk::Reference &chroms,
                         std::uint32_t bin_size, bool overwrite_if_exists,
                         hictk::StandardAttributes attributes) {
  *this = File::create_new_cooler<PixelT>(uri, chroms, bin_size, overwrite_if_exists, attributes);
}

inline void File::close() {
  this->finalize();
  *this = File{};
}

inline void File::finalize() {
  if (!_fp || !_finalize) {
    assert(!_bins == !_fp);
    assert(!_index == !_fp);
    return;
  }

  assert(this->_bins);
  assert(this->_index);
  try {
    this->write_chromosomes();
    this->write_bin_table();

    assert(_attrs.nnz.has_value());
    _index->nnz() = static_cast<std::uint64_t>(*_attrs.nnz);
    this->write_indexes();
    this->write_attributes();

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while closing file {}: {}\n"
                               "File is likely corrupted or incomplete"),
                    this->path(), e.what()));
  }
}

inline HighFive::File File::open_file(std::string_view uri, unsigned int mode, bool validate) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  const auto [file_path, root_grp] = parse_cooler_uri(uri);

  const auto new_file = !std::filesystem::exists(file_path);
  HighFive::File f(file_path, mode);
  if (!validate || new_file) {
    return f;
  }

  const auto status = utils::is_cooler(f, root_grp);
  if (!status) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("\"{}\" does not look like a valid Cooler file:\n"
                               "Validation report:\n{}"),
                    uri, status));
  }

  return f;
}

inline auto File::open_or_create_root_group(HighFive::File &f, std::string_view uri) -> RootGroup {
  if (f.exist(parse_cooler_uri(uri).group_path)) {
    return open_root_group(f, uri);
  }
  return create_root_group(f, uri);
}

namespace internal {
template <typename Variant, std::size_t i = 0>
[[nodiscard]] inline Variant read_pixel_variant(const HighFive::DataSet &dset) {
  if constexpr (i < std::variant_size_v<Variant>) {
    using T = std::variant_alternative_t<i, Variant>;
    if (dset.getDataType() != HighFive::create_datatype<T>()) {
      return read_pixel_variant<Variant, i + 1>(dset);
    }
    return T{};
  }

  constexpr bool variant_has_monostate =
      std::is_same_v<std::monostate, std::variant_alternative_t<0, Variant>>;
  if constexpr (variant_has_monostate) {
    return std::monostate();
  } else {
    throw std::runtime_error(
        fmt::format(FMT_STRING("Unsupported type for dataset \"{}\""), dset.getPath()));
  }
}
}  // namespace internal

inline internal::NumericVariant File::detect_pixel_type(const RootGroup &root_grp,
                                                        std::string_view path) {
  [[maybe_unused]] HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto dset = root_grp().getDataSet(std::string{path});
  return internal::read_pixel_variant<internal::NumericVariant>(dset);
}

template <typename N, bool cis>
inline void File::update_pixel_sum(N partial_sum) {
  static_assert(std::is_arithmetic_v<N>);

  auto &buff = cis ? this->_attrs.cis : this->_attrs.sum;
  assert(buff.has_value());
  if constexpr (std::is_floating_point_v<N>) {
    std::get<double>(*buff) += conditional_static_cast<double>(partial_sum);
  } else {
    assert(std::is_integral_v<N>);
    std::get<std::int64_t>(*buff) += conditional_static_cast<std::int64_t>(partial_sum);
  }
}

}  // namespace hictk
