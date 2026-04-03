// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

// clang-format off
#include "hictk/cooler/cooler.hpp"
// clang-format on

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>
#include <highfive/H5Utility.hpp>
#include <iterator>
#include <string>
#include <string_view>
#include <utility>
#include <variant>

#include "hictk/bin.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/common.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/index.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler {

void File::flush() { _root_group().getFile().flush(); }

auto File::create_root_group(HighFive::File &f, std::string_view uri, bool write_sentinel_attr)
    -> RootGroup {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  auto grp = f.createGroup(parse_cooler_uri(uri).group_path);
  if (write_sentinel_attr) {
    Attribute::write(grp, internal::SENTINEL_ATTR_NAME, internal::SENTINEL_ATTR_VALUE);
    f.flush();
  }

  return RootGroup{grp};
}

auto File::create_groups(RootGroup &root_grp) -> GroupMap {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
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

auto File::create_groups(RootGroup &root_grp, Group chroms_grp, Group bins_grp) -> GroupMap {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
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

void File::write_standard_attributes(RootGroup &root_grp, const Attributes &attributes,
                                     bool skip_sentinel_attr) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  if (attributes.assembly) {
    Attribute::write(root_grp(), "assembly", *attributes.assembly);
  }
  if (attributes.bin_size == 0) {
    assert(attributes.bin_type == BinTable::Type::variable);
    Attribute::write(root_grp(), "bin-size", std::string{"null"});
  } else {
    assert(attributes.bin_type == BinTable::Type::fixed);
    Attribute::write(root_grp(), "bin-size", attributes.bin_size);
  }
  const std::string bin_type = attributes.bin_type == BinTable::Type::fixed ? "fixed" : "variable";
  Attribute::write(root_grp(), "bin-type", bin_type);
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

void File::write_attributes(bool skip_sentinel_attr) {
  assert(_attrs.nbins == bins().size());
  assert(_attrs.nchroms == chromosomes().size());
  assert(_attrs.nnz == _datasets.at("pixels/count").size());

  write_standard_attributes(_root_group, _attrs, skip_sentinel_attr);
  flush();
  if (skip_sentinel_attr) {
    using T [[maybe_unused]] = remove_cvref_t<decltype(internal::SENTINEL_ATTR_VALUE)>;
    assert(Attribute::read<T>(_root_group(), internal::SENTINEL_ATTR_NAME) ==
           internal::SENTINEL_ATTR_VALUE);
    Attribute::write(_root_group(), "format-version", _attrs.format_version, true);
    flush();
  }
}

void File::write_chromosomes() {
  assert(_datasets.contains("chroms/name"));
  assert(_datasets.contains("chroms/length"));
  assert(!chromosomes().empty());

  write_chromosomes(dataset("chroms/name"), dataset("chroms/length"), chromosomes().begin(),
                    chromosomes().end());

  _attrs.nchroms = static_cast<std::int32_t>(chromosomes().size());
}

void File::write_bin_table() {
  write_bin_table(dataset("bins/chrom"), dataset("bins/start"), dataset("bins/end"), bins());

  _attrs.nbins = bins().size();
}

void File::write_bin_table(Dataset &chrom_dset, Dataset &start_dset, Dataset &end_dset,
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

void File::write_indexes() {
  assert(_attrs.nnz.has_value());
  index().finalize(static_cast<std::uint64_t>(*_attrs.nnz));
  write_indexes(dataset("indexes/chrom_offset"), dataset("indexes/bin1_offset"), index());
}

void File::write_indexes(Dataset &chrom_offset_dset, Dataset &bin_offset_dset, const Index &idx) {
  chrom_offset_dset.write(idx.compute_chrom_offsets(), 0, true);

  bin_offset_dset.write(idx.begin(), idx.end(), 0, true);

  assert(chrom_offset_dset.size() == idx.chromosomes().size() + 1);
  assert(bin_offset_dset.size() == idx.size() + 1);
}

void File::write_sentinel_attr(HighFive::Group grp) {
  assert(!check_sentinel_attr(grp));

  Attribute::write(grp, internal::SENTINEL_ATTR_NAME, internal::SENTINEL_ATTR_VALUE, true);
  grp.getFile().flush();
}

void File::write_sentinel_attr() { write_sentinel_attr(_root_group()); }

}  // namespace hictk::cooler
