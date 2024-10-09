// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/btree.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

inline SingleCellAttributes SingleCellAttributes::init(std::uint32_t bin_size_) {
  SingleCellAttributes attrs{};
  attrs.bin_size = bin_size_;
  attrs.bin_type = bin_size_ == 0 ? BinTable::Type::variable : BinTable::Type::fixed;
  return attrs;
}

inline SingleCellAttributes SingleCellAttributes::init_empty() noexcept {
  SingleCellAttributes attrs{};

  attrs.creation_date.reset();
  attrs.format_url.reset();
  attrs.generated_by.reset();
  attrs.assembly.reset();
  attrs.nbins.reset();
  attrs.nchroms.reset();
  attrs.metadata.reset();
  attrs.storage_mode.reset();

  return attrs;
}

inline bool SingleCellAttributes::operator==(const SingleCellAttributes& other) const noexcept {
  // clang-format off
  return bin_size == other.bin_size &&
         bin_type == other.bin_type &&
         format == other.format &&
         format_version == other.format_version &&
         nbins == other.nbins &&
         nchroms == other.nchroms &&
         ncells == other.ncells;
  // clang-format on
}

inline bool SingleCellAttributes::operator!=(const SingleCellAttributes& other) const noexcept {
  return !(*this == other);
}

inline SingleCellFile::SingleCellFile(const HighFive::File& fp, BinTable bins,
                                      SingleCellAttributes attrs)
    : _root_grp(std::make_unique<RootGroup>(RootGroup{fp.getGroup("/")})),
      _cells(read_cells(fp)),
      _attrs(std::move(attrs)),
      _bins(std::make_shared<const BinTable>(std::move(bins))) {
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

inline SingleCellFile::SingleCellFile(const std::filesystem::path& path, unsigned int mode)
    : SingleCellFile(HighFive::File(path.string(), mode),
                     init_bin_table(HighFive::File(path.string())),
                     read_standard_attributes(HighFive::File(path.string()),
                                              mode != HighFive::File::ReadOnly)) {}

inline SingleCellFile SingleCellFile::create(const std::filesystem::path& path,
                                             const Reference& chroms, std::uint32_t bin_size,
                                             bool force_overwrite,
                                             SingleCellAttributes attributes) {
  return SingleCellFile::create(path, BinTable{chroms, bin_size}, force_overwrite,
                                std::move(attributes));
}

inline SingleCellFile SingleCellFile::create(const std::filesystem::path& path, BinTable bins,
                                             bool force_overwrite,
                                             SingleCellAttributes attributes) {
  if (!force_overwrite && std::filesystem::exists(path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to initialize file \"{}\": file already exists"), path));
  }

  if (force_overwrite) {
    std::filesystem::remove(path);
  }

  const HighFive::File fp(path.string(), HighFive::File::Create);
  RootGroup root_grp{fp.getGroup("/")};

  attributes.bin_size = bins.resolution();
  attributes.bin_type = bins.resolution() == 0 ? BinTable::Type::variable : BinTable::Type::fixed;
  create_groups(root_grp);
  create_datasets(root_grp, bins);

  Dataset chrom_name_dset(root_grp, root_grp().getDataSet("chroms/name"));
  Dataset chrom_size_dset(root_grp, root_grp().getDataSet("chroms/length"));
  File::write_chromosomes(chrom_name_dset, chrom_size_dset, bins.chromosomes().begin(),
                          bins.chromosomes().end());
  attributes.nchroms = bins.chromosomes().size();

  Dataset bins_chrom_dset(root_grp, root_grp().getDataSet("bins/chrom"));
  Dataset bins_start_dset(root_grp, root_grp().getDataSet("bins/start"));
  Dataset bins_end_dset(root_grp, root_grp().getDataSet("bins/end"));
  File::write_bin_table(bins_chrom_dset, bins_start_dset, bins_end_dset, bins);
  attributes.nbins = bins.size();

  write_standard_attributes(root_grp, attributes);
  return {fp, std::move(bins), attributes};
}

constexpr const phmap::btree_set<std::string>& SingleCellFile::cells() const noexcept {
  return _cells;
}
constexpr const SingleCellAttributes& SingleCellFile::attributes() const noexcept { return _attrs; }

inline File SingleCellFile::open(std::string_view cell) const {
  assert(_root_grp);
  if (!_cells.contains(cell)) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("unable to find cell \"{}\" in file {}"), cell, path()));
  }
  return File(RootGroup{(*_root_grp)().getGroup(fmt::format(FMT_STRING("cells/{}"), cell))});
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

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

inline SingleCellFile::operator bool() const noexcept { return !!_root_grp; }

inline std::string SingleCellFile::path() const { return (*_root_grp)().getFile().getName(); }

inline auto SingleCellFile::bins() const noexcept -> const BinTable& { return *_bins; }
inline std::uint32_t SingleCellFile::resolution() const noexcept { return bins().resolution(); }
inline auto SingleCellFile::chromosomes() const noexcept -> const Reference& {
  return bins().chromosomes();
}

inline HighFive::File SingleCellFile::file_handle() {
  assert(_root_grp);
  return (*_root_grp)().getFile();
}
inline const HighFive::File& SingleCellFile::file_handle() const {
  assert(_root_grp);
  return (*_root_grp)().getFile();
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

HICTK_DISABLE_WARNING_PUSH HICTK_DISABLE_WARNING_UNREACHABLE_CODE inline SingleCellAttributes
SingleCellFile::read_standard_attributes(const HighFive::File& f, bool initialize_missing) {
  const RootGroup root_grp{f.getGroup("/")};
  auto attrs =
      initialize_missing ? SingleCellAttributes::init(0) : SingleCellAttributes::init_empty();
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT

  auto read_or_throw = [&](const auto& key, auto& buff) {
    using T = remove_cvref_t<decltype(buff)>;
    try {
      buff = Attribute::read<T>(root_grp(), key);
    } catch (const std::exception& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"),
                      key, root_grp().getPath(), e.what()));
    }
  };

  // Read mandatory attributes
  read_or_throw("bin-size", attrs.bin_size);
  std::string bin_type{"fixed"};
  read_or_throw("bin-type", bin_type);
  attrs.bin_type = bin_type == "fixed" ? BinTable::Type::fixed : BinTable::Type::variable;
  read_or_throw("format", attrs.format);
  read_or_throw("format-version", attrs.format_version);

  // Try to read reserved attributes
  auto missing_ok = true;
  internal::read_optional(root_grp, "creation-date", attrs.creation_date, missing_ok);
  internal::read_optional(root_grp, "format-url", attrs.format_url, missing_ok);
  internal::read_optional(root_grp, "generated-by", attrs.generated_by, missing_ok);

  if (!internal::read_optional(root_grp, "genome-assembly", attrs.assembly, missing_ok)) {
    internal::read_optional(root_grp, "assembly", attrs.assembly, missing_ok);
  }

  internal::read_optional(root_grp, "metadata", attrs.metadata, missing_ok);
  internal::read_optional(root_grp, "storage-mode", attrs.storage_mode, missing_ok);

  // Try to read other common attributes
  internal::read_optional(root_grp, "nbins", attrs.nbins, missing_ok);
  internal::read_optional(root_grp, "ncells", attrs.ncells, missing_ok);
  internal::read_optional(root_grp, "nchroms", attrs.nchroms, missing_ok);

  return attrs;
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)
HICTK_DISABLE_WARNING_POP

inline BinTable SingleCellFile::init_bin_table(const HighFive::File& f) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  const RootGroup root_grp{f.getGroup("/")};
  auto chroms = File::import_chroms(Dataset(root_grp, f.getDataSet("/chroms/name")),
                                    Dataset(root_grp, f.getDataSet("/chroms/length")), false);
  const auto bin_type = Attribute::read<std::string>(root_grp(), "bin-type");
  if (bin_type == "fixed") {
    const auto bin_size = Attribute::read<std::uint32_t>(root_grp(), "bin-size");

    return {std::move(chroms), bin_size};
  }

  assert(bin_type == "variable");
  return {std::move(chroms),
          Dataset(root_grp, f.getDataSet("bins/start")).read_all<std::vector<std::uint32_t>>(),
          Dataset(root_grp, f.getDataSet("bins/end")).read_all<std::vector<std::uint32_t>>()};
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

inline phmap::btree_set<std::string> SingleCellFile::read_cells(const HighFive::File& f) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  phmap::btree_set<std::string> cells{};
  auto buffer = f.getGroup("/cells").listObjectNames();
  std::move(buffer.begin(), buffer.end(), std::inserter(cells, cells.begin()));
  return cells;
}

inline void SingleCellFile::create_groups(RootGroup& root_grp) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto bins_group = root_grp().createGroup("/bins");
  auto chroms_group = root_grp().createGroup("/chroms");
  auto res_group = root_grp().createGroup("/cells");
}

inline void SingleCellFile::write_standard_attributes(RootGroup& root_grp,
                                                      const SingleCellAttributes& attrs) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  Attribute::write(root_grp(), "assembly",
                   attrs.assembly.has_value() ? *attrs.assembly : "unknown");
  Attribute::write(root_grp(), "bin-size", attrs.bin_size);
  const std::string bin_type = attrs.bin_type == BinTable::Type::fixed ? "fixed" : "variable";
  Attribute::write(root_grp(), "bin-type", bin_type);
  if (attrs.creation_date.has_value()) {
    Attribute::write(root_grp(), "creation-date", *attrs.creation_date);
  }
  Attribute::write(root_grp(), "format", attrs.format);
  if (attrs.format_url.has_value()) {
    Attribute::write(root_grp(), "format-url", *attrs.format_url);
  }
  Attribute::write(root_grp(), "format-version", static_cast<std::int64_t>(attrs.format_version));
  if (attrs.generated_by.has_value()) {
    Attribute::write(root_grp(), "generated-by", *attrs.generated_by);
  }
  if (attrs.metadata.has_value()) {
    Attribute::write(root_grp(), "metadata", *attrs.metadata);
  }
  assert(attrs.nbins.has_value());
  assert(attrs.nchroms.has_value());
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  Attribute::write(root_grp(), "nbins", *attrs.nbins);
  Attribute::write(root_grp(), "ncells", 0);
  // NOLINTNEXTLINE(bugprone-unchecked-optional-access)
  Attribute::write(root_grp(), "nchroms", *attrs.nchroms);
  if (attrs.storage_mode.has_value()) {
    Attribute::write(root_grp(), "storage-mode", *attrs.storage_mode);
  }
}

inline void SingleCellFile::create_datasets(RootGroup& root_grp, const BinTable& bins) {
  const auto default_aprop =
      Dataset::init_access_props(DEFAULT_HDF5_CHUNK_SIZE, DEFAULT_HDF5_DATASET_CACHE_SIZE, 1.0);

  auto create_dataset = [&](const auto& path, const auto& type, auto aprop) {
    using T = remove_cvref_t<decltype(type)>;
    if constexpr (is_string_v<T>) {
      const auto& chrom_with_longest_name = bins.chromosomes().chromosome_with_longest_name();
      Dataset{root_grp, path, chrom_with_longest_name.name(), HighFive::DataSpace::UNLIMITED,
              aprop};
    } else {
      Dataset{root_grp, path, type, HighFive::DataSpace::UNLIMITED, aprop};
    }
  };

  create_dataset("chroms/name", std::string{}, default_aprop);
  create_dataset("chroms/length", std::int32_t{}, default_aprop);

  create_dataset("bins/chrom", std::int32_t{}, default_aprop);
  create_dataset("bins/start", std::int32_t{}, default_aprop);
  create_dataset("bins/end", std::int32_t{}, default_aprop);
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

  auto create_dataset = [&](const auto& path, const auto& type, auto aprop, auto cprop) {
    Dataset{root_grp, path, type, HighFive::DataSpace::UNLIMITED, aprop, cprop};
  };

  create_dataset("pixels/bin1_id", std::int64_t{}, pixels_aprop, default_cprop);
  create_dataset("pixels/bin2_id", std::int64_t{}, pixels_aprop, default_cprop);
  create_dataset("pixels/count", PixelT{}, pixels_aprop, default_cprop);

  create_dataset("indexes/bin1_offset", std::int64_t{}, default_aprop, default_cprop);
  create_dataset("indexes/chrom_offset", std::int64_t{}, default_aprop, default_cprop);
}

}  // namespace hictk::cooler
