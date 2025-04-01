// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/chrono.h>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <fmt/std.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <highfive/H5Exception.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Utility.hpp>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/common.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/common.hpp"
#include "hictk/cooler/dataset.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/uri.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::cooler::utils {

constexpr ValidationStatusCooler::operator bool() const noexcept { return is_cooler; }

constexpr ValidationStatusMultiresCooler::operator bool() const noexcept {
  return is_multires_file;
}

constexpr ValidationStatusScool::operator bool() const noexcept { return is_scool_file; }

inline ValidationStatusCooler is_cooler(std::string_view uri) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  const auto [file_path, root_path] = parse_cooler_uri(uri);
  try {
    const HighFive::File fp(file_path, HighFive::File::ReadOnly);
    return is_cooler(fp, root_path);
  } catch (const std::exception &e) {
    const std::string_view msg{e.what()};
    ValidationStatusCooler s{};
    s.uri = std::string{uri};
    s.unable_to_open_file = msg.find("Unable to open file") != std::string_view::npos;
    s.is_hdf5 = msg.find("Not an HDF5 file") == std::string_view::npos;

    if (!s.unable_to_open_file && !s.is_hdf5) {
      throw;
    }
    return s;
  }
}

inline ValidationStatusMultiresCooler is_multires_file(std::string_view uri,
                                                       bool validate_resolutions, bool exhaustive,
                                                       std::int64_t min_version) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  try {
    const HighFive::File fp(std::string{uri}, HighFive::File::ReadOnly);
    return is_multires_file(fp, validate_resolutions, exhaustive, min_version);
  } catch (const std::exception &e) {
    const std::string_view msg{e.what()};
    ValidationStatusMultiresCooler s{};
    s.uri = std::string{uri};
    s.unable_to_open_file = msg.find("Unable to open file") != std::string_view::npos;
    s.is_hdf5 = msg.find("Not an HDF5 file") == std::string_view::npos;

    if (!s.unable_to_open_file && !s.is_hdf5) {
      throw;
    }
    return s;
  }
}

inline ValidationStatusScool is_scool_file(std::string_view uri, bool validate_cells,
                                           bool exhaustive) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  try {
    const HighFive::File fp(std::string{uri}, HighFive::File::ReadOnly);
    return is_scool_file(fp, validate_cells, exhaustive);
  } catch (const std::exception &e) {
    const std::string_view msg{e.what()};
    ValidationStatusScool s{};
    s.uri = std::string{uri};
    s.unable_to_open_file = msg.find("Unable to open file") != std::string_view::npos;
    s.is_hdf5 = msg.find("Not an HDF5 file") == std::string_view::npos;

    if (!s.unable_to_open_file && !s.is_hdf5) {
      throw;
    }
    return s;
  }
}

inline ValidationStatusCooler is_cooler(const HighFive::File &fp, std::string_view root_path) {
  try {
    return is_cooler(fp.getGroup(std::string{root_path}));
  } catch (const HighFive::Exception &) {
    [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT

    const auto root_group = fp.getGroup("/");

    ValidationStatusCooler status{};
    status.unable_to_open_file = false;
    status.is_hdf5 = fp.isValid();
    if (Attribute::exists(root_group, "format-version")) {
      const auto version = Attribute::read<std::uint8_t>(root_group, "format-version");
      status.file_was_properly_closed = version != cooler::internal::SENTINEL_ATTR_VALUE;
    }

    if (!root_path.empty() && root_path != "/") {
      // NOLINTBEGIN(readability-implicit-bool-conversion)
      status.uri += fmt::format(FMT_STRING("{}::/{}"), fp.getName(),
                                root_path.substr(root_path.front() == '/'));
      status.missing_groups.emplace_back(root_path.substr(root_path.front() == '/'));
      // NOLINTEND(readability-implicit-bool-conversion)
    }

    status.is_cooler = false;

    return status;
  }
}

inline ValidationStatusCooler is_cooler(const HighFive::Group &root_group) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  ValidationStatusCooler status{};
  status.uri = root_group.getFile().getName();
  const auto root_path = root_group.getPath();
  if (!root_path.empty() && root_path != "/") {
    // NOLINTNEXTLINE(readability-implicit-bool-conversion)
    status.uri += fmt::format(FMT_STRING("::/{}"), root_path.substr(root_path.front() == '/'));
  }

  status.is_hdf5 = root_group.getFile().isValid();
  // Check file is in HDF5 format
  if (!status.is_hdf5) {
    return status;
  }

  // Check file has the appropriate format attribute
  if (Attribute::exists(root_group, "format")) {
    const auto format = Attribute::read<std::string>(root_group, "format");
    status.missing_or_invalid_format_attr = format != COOL_MAGIC;
  }

  if (Attribute::exists(root_group, "format-version")) {
    const auto version = Attribute::read<std::uint8_t>(root_group, "format-version");
    status.file_was_properly_closed = version != cooler::internal::SENTINEL_ATTR_VALUE;
    status.missing_or_invalid_format_attr |= version == 0 || version > 3;
  }

  // Check file has a bin-type that we support
  if (Attribute::exists(root_group, "bin-type")) {
    const auto bin_type = Attribute::read<std::string>(root_group, "bin-type");
    status.missing_or_invalid_bin_type_attr = bin_type != "fixed" && bin_type != "variable";
  }

  // Check file has the mandatory groups
  for (const auto &group : MANDATORY_GROUP_NAMES) {
    if (!root_group.exist(std::string{group})) {
      status.missing_groups.emplace_back(group);
      continue;
    }
    try {
      root_group.getGroup(std::string{group});
    } catch (const HighFive::Exception &) {
      status.missing_groups.emplace_back(group);
    }
  }

  // clang-format off
  status.is_cooler = status.is_hdf5 &&
                     status.file_was_properly_closed &&
                     !status.missing_or_invalid_format_attr &&
                     !status.missing_or_invalid_bin_type_attr &&
                     status.missing_groups.empty();
  // clang-format on

  return status;
}

inline ValidationStatusMultiresCooler is_multires_file(const HighFive::File &fp,
                                                       bool validate_resolutions, bool exhaustive,
                                                       std::int64_t min_version) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  ValidationStatusMultiresCooler status{};
  status.uri = fp.getName();

  status.file_was_properly_closed = true;
  status.is_hdf5 = fp.isValid();
  // Check file is in HDF5 format
  if (!status.is_hdf5) {
    return status;
  }

  // Check file has the appropriate format attribute
  if (Attribute::exists(fp, "format")) {
    const auto format = Attribute::read<std::string>(fp, "format");
    status.missing_or_invalid_format_attr = format != MCOOL_MAGIC;
  }

  // Check file has the appropriate format attribute
  if (Attribute::exists(fp, "format-version")) {
    const auto version = Attribute::read<std::uint8_t>(fp, "format-version");
    status.missing_or_invalid_format_attr |= version == 0 || version > 3;
  }

  // Check file has a bin-type that we support
  // NOTE: .mcool files are not required to advertise the bin type they are using at the root level
  status.missing_or_invalid_bin_type_attr = false;
  if (Attribute::exists(fp, "bin-type")) {
    const auto bin_type = Attribute::read<std::string>(fp, "bin-type");
    status.missing_or_invalid_bin_type_attr = bin_type != "fixed" && bin_type != "variable";
  }

  // Try to read resolutions from the Cooler's root
  const auto resolutions = [&]() -> std::vector<std::string> {
    try {
      auto root_grp = fp.getGroup("/resolutions");
      auto resolutions_ = root_grp.listObjectNames();

      // TODO: I am not sure I fully understand the docs. here:
      //       https://cooler.readthedocs.io/en/latest/schema.html#multi-resolution
      if (min_version >= 2) {
        return resolutions_;
      }

      // The old multi-resolution layout used resolutions strictly in increments of powers of 2.
      // In this layout (MCOOL version 2), the data collections are named by zoom level,
      // starting with XYZ.1000.mcool::0 being the coarsest resolution up until the finest or
      // “base” resolution (e.g., XYZ.1000.mcool::14 for 14 levels of coarsening).
      auto it = std::find(resolutions_.begin(), resolutions_.end(), "0");
      if (it != resolutions_.end()) {
        return {};
      }
      return resolutions_;
    } catch (const HighFive::GroupException &) {
      return {};
    }
  }();

  if (resolutions.empty()) {
    status.missing_groups.emplace_back("resolutions");
  }

  if (validate_resolutions) {
    for (const auto &resolution : resolutions) {
      const auto suffix = fmt::format(FMT_STRING("resolutions/{}"), resolution);

      if (auto status_ = is_cooler(fp, suffix); !status_) {
        status.file_was_properly_closed &= status_.file_was_properly_closed;
        status.invalid_resolutions.emplace(resolution, std::move(status_));
        if (!exhaustive) {
          break;
        }
      } else if (exhaustive) {
        status.valid_resolutions.emplace(
            hictk::internal::parse_numeric_or_throw<std::uint32_t>(resolution), std::move(status_));
      }
    }
  }

  // clang-format off
  status.is_multires_file = status.is_hdf5 &&
                            status.file_was_properly_closed &&
                            !status.missing_or_invalid_format_attr &&
                            !status.missing_or_invalid_bin_type_attr &&
                            status.missing_groups.empty() &&
                            status.invalid_resolutions.empty();
  // clang-format on

  return status;
}

inline ValidationStatusScool is_scool_file(const HighFive::File &fp, bool validate_cells,
                                           bool exhaustive) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  ValidationStatusScool status{};
  status.uri = fp.getName();

  status.file_was_properly_closed = true;
  status.is_hdf5 = fp.isValid();
  // Check file is in HDF5 format
  if (!status.is_hdf5) {
    return status;
  }

  // Check file has the appropriate format attribute
  if (Attribute::exists(fp, "format")) {
    const auto format = Attribute::read<std::string>(fp, "format");
    status.missing_or_invalid_format_attr = format != SCOOL_MAGIC;
  }

  if (Attribute::exists(fp, "format")) {
    const auto version = Attribute::read<std::uint8_t>(fp, "format-version");
    status.missing_or_invalid_format_attr |= version == 0 || version > 3;
  }

  // Check file has a bin-type that we support
  // NOTE: .scool files are not required to advertise the bin type they are using at the root level
  status.missing_or_invalid_bin_type_attr = false;
  if (Attribute::exists(fp, "bin-type")) {
    const auto bin_type = Attribute::read<std::string>(fp, "bin-type");
    status.missing_or_invalid_bin_type_attr = bin_type != "fixed" && bin_type != "variable";
  }

  constexpr std::array<std::string_view, 3> scool_root_groups{"chroms", "bins", "cells"};

  const auto root_groups = fp.listObjectNames();

  // Check file has the mandatory groups
  for (const auto &group : scool_root_groups) {
    try {
      fp.getGroup(std::string{group});
    } catch (const HighFive::Exception &) {
      status.missing_groups.emplace_back(group);
    }
  }

  const auto cells = [&]() -> std::vector<std::string> {
    try {
      auto root_grp = fp.getGroup("/cells");
      return root_grp.listObjectNames();
    } catch (const HighFive::GroupException &) {
      return {};
    }
  }();

  // Try to check if the number of groups under /cells is consistent with available metadata
  status.unexpected_number_of_cells = false;
  if (Attribute::exists(fp, "ncells")) {
    const auto expected_num_cells = Attribute::read<std::uint64_t>(fp, "ncells");
    status.unexpected_number_of_cells =
        expected_num_cells != conditional_static_cast<std::uint64_t>(cells.size());
  }

  if (validate_cells) {
    for (const auto &cell : cells) {
      const auto suffix = fmt::format(FMT_STRING("cells/{}"), cell);
      if (auto status_ = is_cooler(fp, suffix); !status_) {
        status.file_was_properly_closed &= status_.file_was_properly_closed;
        status.invalid_cells.emplace(cell, std::move(status_));
        if (!exhaustive) {
          break;
        }
      } else if (exhaustive) {
        status.valid_cells.emplace(cell, std::move(status_));
      }
    }
  }

  // clang-format off
  status.is_scool_file = status.is_hdf5 &&
                         status.file_was_properly_closed &&
                         !status.missing_or_invalid_format_attr &&
                         !status.missing_or_invalid_bin_type_attr &&
                         status.missing_groups.empty() &&
                         !status.unexpected_number_of_cells &&
                         status.invalid_cells.empty();
  // clang-format on

  return status;
}

[[nodiscard]] inline bool index_is_valid(std::string_view uri) {
  [[maybe_unused]] std::string buff{};
  return index_is_valid(uri, buff);
}

[[nodiscard]] inline bool index_is_valid(std::string_view uri, std::string &error_buffer) {
  // See https://github.com/robomics/20221129_4dnucleome_bug_report
  // and https://github.com/open2c/cooler/issues/319

  error_buffer.clear();
  if (!is_cooler(uri)) {
    return false;
  }

  const auto clr_uri = parse_cooler_uri(uri);

  const HighFive::File h5(clr_uri.file_path);

  const auto root_grp = h5.getGroup(clr_uri.group_path.empty() ? "/" : clr_uri.group_path);

  // check storage-mode attribute is valid
  std::visit(
      [&](const auto &storage_mode) {
        using T = remove_cvref_t<decltype(storage_mode)>;
        if constexpr (std::is_same_v<T, std::monostate>) {
          return;
        }
        if constexpr (is_string_v<T>) {
          if (storage_mode != "symmetric-upper") {
            throw std::runtime_error(fmt::format(
                FMT_STRING(
                    "validating the index of Coolers with storage-mode=\"{}\" is not supported"),
                storage_mode));
          }

          return;
        }

        throw std::runtime_error(
            fmt::format(FMT_STRING("unknown storage-mode=\"{}\""), storage_mode));
      },
      Attribute::read(root_grp, "storage-mode", true));

  if (!root_grp.exist("indexes/bin1_offset")) {
    throw std::runtime_error("file does not have dataset \"indexes/bin1_offset\"");
  }
  if (!root_grp.exist("pixels/bin2_id")) {
    throw std::runtime_error("file does not have dataset \"pixels/bin2_id\"");
  }

  const cooler::Dataset bin1_dset{RootGroup{root_grp}, root_grp.getDataSet("indexes/bin1_offset")};
  const cooler::Dataset bin2_dset{RootGroup{root_grp}, root_grp.getDataSet("pixels/bin2_id")};

  auto is_sorted = [&](auto &first, const auto &last) {
    // we use this instead of std::is_sorted for performance reasons.
    // Basically, we want to pass iterators by reference such that when advancing first, we never
    // need to read a chunk more than once
    assert(first < last);

    auto x1 = *first;
    std::ignore = ++first;
    while (first != last) {
      const auto x2 = *first;
      if (x1 > x2) {
        return false;
      }
      x1 = x2;
      std::ignore = ++first;
    }

    return true;
  };

  auto first = bin2_dset.begin<std::uint64_t>();
  const auto bin1_offset = bin1_dset.read_all<std::vector<std::uint64_t>>();
  assert(!bin1_offset.empty());
  for (std::size_t i1 = 1; i1 < bin1_offset.size() - 1; ++i1) {
    const auto i0 = i1 - 1;

    const auto j0 = bin1_offset[i0];
    const auto j1 = bin1_offset[i1];

    if (j0 + 1 >= j1) {
      // no point in checking empty rows or rows with only 1 pixel
      continue;
    }

    first.seek(j0);
    const auto last = bin2_dset.make_iterator_at_offset<std::uint64_t>(j1, 0);

    if (!is_sorted(first, last)) {
      error_buffer = fmt::format(FMT_STRING("pixels between {}-{} (bin1_id={}) are not sorted in "
                                            "ascending order (and very likely "
                                            "contain duplicate entries)"),
                                 bin1_offset[i0], bin1_offset[i1], i0);
      return false;
    }
  }

  return true;
}

}  // namespace hictk::cooler::utils

constexpr auto fmt::formatter<hictk::cooler::utils::ValidationStatusCooler>::parse(
    format_parse_context &ctx) -> format_parse_context::iterator {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw format_error("invalid format");
  }
  return ctx.end();
}

inline auto fmt::formatter<hictk::cooler::utils::ValidationStatusCooler>::format(
    const hictk::cooler::utils::ValidationStatusCooler &s, format_context &ctx)
    -> decltype(ctx.out()) {
  // clang-format off
  return fmt::format_to(
      ctx.out(),
      FMT_STRING("uri=\"{}\"\n"
                 "is_hdf5={}\n"
                 "unable_to_open_file={}\n"
                 "file_was_properly_closed={}\n"
                 "missing_or_invalid_format_attr={}\n"
                 "missing_or_invalid_bin_type_attr={}\n"
                 "missing_groups=[{}]\n"
                 "is_valid_cooler={}"),
      s.uri,
      s.is_hdf5,
      s.unable_to_open_file,
      s.file_was_properly_closed,
      s.missing_or_invalid_format_attr,
      s.missing_or_invalid_bin_type_attr,
      fmt::join(s.missing_groups, ", "),
      s.is_cooler);
  // clang-format on
}

constexpr auto fmt::formatter<hictk::cooler::utils::ValidationStatusMultiresCooler>::parse(
    format_parse_context &ctx) -> format_parse_context::iterator {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw format_error("invalid format");
  }
  return ctx.end();
}

inline auto fmt::formatter<hictk::cooler::utils::ValidationStatusMultiresCooler>::format(
    const hictk::cooler::utils::ValidationStatusMultiresCooler &s, format_context &ctx)
    -> decltype(ctx.out()) {
  // clang-format off
  return fmt::format_to(
      ctx.out(),
      FMT_STRING("uri=\"{}\"\n"
                 "is_hdf5={}\n"
                 "unable_to_open_file={}\n"
                 "file_was_properly_closed={}\n"
                 "missing_or_invalid_format_attr={}\n"
                 "missing_or_invalid_bin_type_attr={}\n"
                 "missing_groups=[{}]\n"
                 "is_valid_multires_file={}\n"
                 "invalid_resolutions{}{}"),
      s.uri,
      s.is_hdf5,
      s.unable_to_open_file,
      s.file_was_properly_closed,
      s.missing_or_invalid_format_attr,
      s.missing_or_invalid_bin_type_attr,
      fmt::join(s.missing_groups, ", "),
      s.is_multires_file,
      s.invalid_resolutions.empty() ? "=[]" : ":\n - ",
      fmt::join(s.invalid_resolutions, "\n - "));
  // clang-format on
}

constexpr auto fmt::formatter<hictk::cooler::utils::ValidationStatusScool>::parse(
    format_parse_context &ctx) -> format_parse_context::iterator {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw format_error("invalid format");
  }
  return ctx.end();
}

inline auto fmt::formatter<hictk::cooler::utils::ValidationStatusScool>::format(
    const hictk::cooler::utils::ValidationStatusScool &s, format_context &ctx)
    -> decltype(ctx.out()) {
  // clang-format off
  return fmt::format_to(
      ctx.out(),
      FMT_STRING("uri=\"{}\"\n"
                 "unable_to_open_file={}\n"
                 "is_hdf5={}\n"
                 "file_was_properly_closed={}\n"
                 "missing_or_invalid_format_attr={}\n"
                 "missing_or_invalid_bin_type_attr={}\n"
                 "missing_groups=[{}]\n"
                 "is_valid_scool_file={}\n"
                 "invalid_cells{}{}"),
      s.uri,
      s.unable_to_open_file,
      s.is_hdf5,
      s.file_was_properly_closed,
      s.missing_or_invalid_format_attr,
      s.missing_or_invalid_bin_type_attr,
      fmt::join(s.missing_groups, ", "),
      s.is_scool_file,
      s.invalid_cells.empty() ? "=[]" : ":\n - ",
      fmt::join(s.invalid_cells, "\n - "));
  // clang-format on
}
