// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <filesystem>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/numeric_variant.hpp"
#include "hictk/reference.hpp"
#include "hictk/suppress_warnings.hpp"

namespace {
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_UNREACHABLE_CODE
template <typename Variant, std::size_t i = 0>
[[nodiscard]] Variant read_pixel_variant(const HighFive::DataSet &dset) {
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
HICTK_DISABLE_WARNING_POP
}  // namespace

namespace hictk::cooler {

File::File(RootGroup entrypoint, HighFiveAccessMode mode, std::size_t cache_size_bytes, double w0,
           bool validate)
    : _mode(mode),
      _root_group(std::move(entrypoint)),
      _groups(open_groups(_root_group)),
      _datasets(open_datasets(_root_group, cache_size_bytes, w0)),
      _attrs(read_standard_attributes(_root_group)),
      _pixel_variant(detect_pixel_type(_root_group)),
      _bins(
          std::make_shared<BinTable>(init_bin_table(_datasets, _attrs.bin_type, _attrs.bin_size))),
      _index(std::make_shared<Index>(init_index(_datasets.at("indexes/chrom_offset"),
                                                _datasets.at("indexes/bin1_offset"), _bins,
                                                _datasets.at("pixels/count").size(), false))) {
  assert(mode == HighFive::File::ReadOnly || mode == HighFive::File::ReadWrite);
  if (validate) {
    validate_bins();
  }
}

File::File(std::string_view uri, std::size_t cache_size_bytes, bool validate)
    : File(open_or_create_root_group(open_file(uri, HighFive::File::ReadOnly, validate), uri),
           HighFive::File::ReadOnly, cache_size_bytes, DEFAULT_HDF5_CACHE_W0, validate) {}

File::File(RootGroup entrypoint, std::size_t cache_size_bytes, bool validate)
    : File(std::move(entrypoint), HighFive::File::ReadOnly, cache_size_bytes, DEFAULT_HDF5_CACHE_W0,
           validate) {}

File File::open_random_access(std::string_view uri, std::size_t cache_size_bytes, bool validate) {
  return File(uri, cache_size_bytes, validate);
}

File File::open_read_once(std::string_view uri, std::size_t cache_size_bytes, bool validate) {
  return {open_or_create_root_group(open_file(uri, HighFive::File::ReadOnly, validate), uri),
          HighFive::File::ReadOnly, cache_size_bytes, 1.0, validate};
}

File File::open_random_access(RootGroup entrypoint, std::size_t cache_size_bytes, bool validate) {
  return File(std::move(entrypoint), cache_size_bytes, validate);
}

File File::open_read_once(RootGroup entrypoint, std::size_t cache_size_bytes, bool validate) {
  return {std::move(entrypoint), HighFive::File::ReadOnly, cache_size_bytes, 1.0, validate};
}

// We need to explicitly define the move ctor because older compilers are not able to automatically
// generate it
File::File(File &&other) noexcept
    : _mode(other._mode),
      _root_group(std::move(other._root_group)),
      _groups(std::move(other._groups)),
      _datasets(std::move(other._datasets)),
      _weights(std::move(other._weights)),
      _weights_scaled(std::move(other._weights_scaled)),
      _attrs(std::move(other._attrs)),
      _pixel_variant(other._pixel_variant),
      _bins(std::move(other._bins)),
      _index(std::move(other._index)),
      _finalize(other._finalize) {
  other._finalize = false;
}

// NOLINTNEXTLINE(bugprone-exception-escape)
File::~File() noexcept {
  try {
    finalize();
  } catch (const std::exception &e) {
    SPDLOG_ERROR(FMT_STRING("an error occurred while closing file \"{}\". File is likely "
                            "corrupted or incomplete. Reason: {}"),
                 path(), e.what());
  } catch (...) {
    SPDLOG_ERROR(FMT_STRING("an unknown error occurred while closing file \"{}\". File is likely "
                            "corrupted or incomplete."),
                 path());
  }
}

// We need to explicitly define the move ctor because older compilers are not able to automatically
// generate it
// NOLINTNEXTLINE(bugprone-exception-escape)
File &File::operator=(File &&other) noexcept {
  if (this == &other) {
    return *this;
  }

  if (_finalize) {
    try {
      finalize();
    } catch (const std::exception &e) {
      SPDLOG_ERROR(FMT_STRING("an error occurred while closing file {}. File is likely "
                              "corrupted or incomplete. Reason: {}"),
                   path(), e.what());
    } catch (...) {
      SPDLOG_ERROR(FMT_STRING("an unknown error occurred while closing file {}. File is likely "
                              "corrupted or incomplete."),
                   path());
    }
  }

  _mode = other._mode;
  _root_group = std::move(other._root_group);
  _groups = std::move(other._groups);
  _datasets = std::move(other._datasets);
  _weights = std::move(other._weights);
  _weights_scaled = std::move(other._weights_scaled);
  _attrs = std::move(other._attrs);
  _pixel_variant = other._pixel_variant;
  _bins = std::move(other._bins);
  _index = std::move(other._index);
  _finalize = other._finalize;

  other._finalize = false;

  return *this;
}

void File::close() {
  finalize();
  _datasets.clear();
  _groups.clear();
  _weights.clear();
  _weights_scaled.clear();
  _pixel_variant = std::int32_t{};
  _bins.reset();
  _index.reset();
  _finalize = false;
  _root_group = RootGroup{};
}

void File::finalize() {
  if (!_bins || !_finalize) {
    assert(!_index == !_bins);
    return;
  }

  assert(_bins);
  assert(_index);
  try {
    assert(_attrs.nnz.has_value());
    _index->set_nnz(static_cast<std::uint64_t>(*_attrs.nnz));
    write_indexes();
    write_attributes();

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("The following error occurred while finalizing file {}: {}\n"
                               "File is likely corrupted or incomplete"),
                    path(), e.what()));
  }
  _finalize = false;
}

HighFive::File File::open_file(std::string_view uri, HighFiveAccessMode mode, bool validate) {
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

auto File::open_or_create_root_group(HighFive::File f, std::string_view uri) -> RootGroup {
  if (f.exist(parse_cooler_uri(uri).group_path)) {
    return open_root_group(f, uri);
  }
  return create_root_group(f, uri);
}

hictk::internal::NumericVariant File::detect_pixel_type(const RootGroup &root_grp,
                                                        std::string_view path) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto dset = root_grp().getDataSet(std::string{path});
  return read_pixel_variant<hictk::internal::NumericVariant>(dset);
}

}  // namespace hictk::cooler
