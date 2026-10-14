// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/singlecell_cooler.hpp"

#include <fmt/format.h>
#include <fmt/std.h>
#include <parallel_hashmap/btree.h>

#include <cassert>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <highfive/H5Utility.hpp>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/pixel_selector.hpp"
#include "hictk/reference.hpp"
#include "hictk/suppress_warnings.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

SingleCellAttributes SingleCellAttributes::init(std::uint32_t bin_size_) {
  SingleCellAttributes attrs{};
  attrs.bin_size = bin_size_;
  attrs.bin_type = bin_size_ == 0 ? BinTable::Type::variable : BinTable::Type::fixed;
  return attrs;
}

SingleCellAttributes SingleCellAttributes::init_empty() noexcept {
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

bool SingleCellAttributes::operator==(const SingleCellAttributes& other) const noexcept {
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

bool SingleCellAttributes::operator!=(const SingleCellAttributes& other) const noexcept {
  return !(*this == other);
}

SingleCellFile::SingleCellFile(const HighFive::File& fp, BinTable bins, SingleCellAttributes attrs)
    : _root_grp(std::make_unique<RootGroup>(RootGroup{fp.getGroup("/")})),
      _cells(read_cells(fp)),
      _attrs(std::move(attrs)),
      _bins(std::make_shared<const BinTable>(std::move(bins))) {
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

SingleCellFile::SingleCellFile(const std::filesystem::path& path, HighFiveAccessMode mode)
    : SingleCellFile(HighFive::File(path.string(), mode),
                     init_bin_table(HighFive::File(path.string())),
                     read_standard_attributes(HighFive::File(path.string()),
                                              mode != HighFive::File::ReadOnly)) {}

SingleCellFile SingleCellFile::create(const std::filesystem::path& path, const Reference& chroms,
                                      std::uint32_t bin_size, bool force_overwrite,
                                      SingleCellAttributes attributes) {
  return SingleCellFile::create(path, BinTable{chroms, bin_size}, force_overwrite,
                                std::move(attributes));
}

SingleCellFile SingleCellFile::create(const std::filesystem::path& path, BinTable bins,
                                      bool force_overwrite, SingleCellAttributes attributes) {
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

File SingleCellFile::open(std::string_view cell) const {
  assert(_root_grp);
  if (!_cells.contains(cell)) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("unable to find cell \"{}\" in file {}"), cell, path()));
  }
  return File(RootGroup{(*_root_grp)().getGroup(fmt::format(FMT_STRING("cells/{}"), cell))});
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

SingleCellFile::operator bool() const noexcept { return !!_root_grp; }

std::string SingleCellFile::path() const { return (*_root_grp)().getFile().getName(); }

auto SingleCellFile::bins() const noexcept -> const BinTable& { return *bins_ptr(); }
auto SingleCellFile::bins_ptr() const noexcept -> std::shared_ptr<const BinTable> { return _bins; }
std::uint32_t SingleCellFile::resolution() const noexcept { return bins().resolution(); }
auto SingleCellFile::chromosomes() const noexcept -> const Reference& {
  return bins().chromosomes();
}

HighFive::File SingleCellFile::file_handle() {
  assert(_root_grp);
  return (*_root_grp)().getFile();
}
const HighFive::File& SingleCellFile::file_handle() const {
  assert(_root_grp);
  return (*_root_grp)().getFile();
}

HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_UNREACHABLE_CODE
SingleCellAttributes SingleCellFile::read_standard_attributes(const HighFive::File& f,
                                                              bool initialize_missing) {
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
  std::ignore = root_grp.try_read_attribute("creation-date", attrs.creation_date);
  std::ignore = root_grp.try_read_attribute("format-url", attrs.format_url);
  std::ignore = root_grp.try_read_attribute("generated-by", attrs.generated_by);

  if (!root_grp.try_read_attribute("genome-assembly", attrs.assembly)) {
    std::ignore = root_grp.try_read_attribute("assembly", attrs.assembly);
  }

  std::ignore = root_grp.try_read_attribute("metadata", attrs.metadata);
  std::ignore = root_grp.try_read_attribute("storage-mode", attrs.storage_mode);

  // Try to read other common attributes
  std::ignore = root_grp.try_read_attribute("nbins", attrs.nbins);
  std::ignore = root_grp.try_read_attribute("ncells", attrs.ncells);
  std::ignore = root_grp.try_read_attribute("nchroms", attrs.nchroms);

  return attrs;
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)
HICTK_DISABLE_WARNING_POP

BinTable SingleCellFile::init_bin_table(const HighFive::File& f) {
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

phmap::btree_set<std::string> SingleCellFile::read_cells(const HighFive::File& f) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  phmap::btree_set<std::string> cells{};
  auto buffer = f.getGroup("/cells").listObjectNames();
  std::move(buffer.begin(), buffer.end(), std::inserter(cells, cells.begin()));
  return cells;
}

void SingleCellFile::create_groups(RootGroup& root_grp) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto bins_group = root_grp().createGroup("/bins");
  auto chroms_group = root_grp().createGroup("/chroms");
  auto res_group = root_grp().createGroup("/cells");
}

void SingleCellFile::write_standard_attributes(RootGroup& root_grp,
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

void SingleCellFile::create_datasets(RootGroup& root_grp, const BinTable& bins) {
  const auto default_aprop =
      Dataset::init_access_props(DEFAULT_HDF5_CHUNK_SIZE, DEFAULT_HDF5_DATASET_CACHE_SIZE, 1.0);

  auto create_dataset = [&](const auto& path, const auto& type, const auto& aprop) {
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

}  // namespace hictk::cooler
