// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/cooler/multires_cooler.hpp"

#include <fmt/format.h>
#include <fmt/std.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <highfive/H5File.hpp>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/group.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/numeric_utils.hpp"
#include "hictk/reference.hpp"

namespace hictk::cooler {

bool MultiResAttributes::operator==(const MultiResAttributes& other) const noexcept {
  return format == other.format && format_version == other.format_version &&
         bin_type == other.bin_type;
}
bool MultiResAttributes::operator!=(const MultiResAttributes& other) const noexcept {
  return !(*this == other);
}

MultiResFile::MultiResFile(const HighFive::File& fp, Reference chroms,
                           std::vector<std::uint32_t> resolutions, MultiResAttributes attrs)
    : _root_grp(std::make_unique<RootGroup>(RootGroup{fp.getGroup("/")})),
      _resolutions(std::move(resolutions)),
      _attrs(std::move(attrs)),
      _chroms(std::move(chroms)) {
  if (_chroms.empty() && !_resolutions.empty()) {
    _chroms = open(_resolutions.back()).chromosomes();
  }
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

MultiResFile::MultiResFile(const std::filesystem::path& path, HighFiveAccessMode mode)
    : MultiResFile(HighFive::File(path.string(), mode), {},
                   read_resolutions(HighFive::File(path.string())),
                   read_attributes(HighFive::File(path.string()))) {}

MultiResFile MultiResFile::create(const std::filesystem::path& path, const Reference& chroms,
                                  bool force_overwrite) {
  if (!force_overwrite && std::filesystem::exists(path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to initialize file \"{}\": file already exists"), path));
  }

  if (force_overwrite) {
    std::filesystem::remove(path);
  }
  HighFive::File fp(path.string(), HighFive::File::Create);

  const MultiResAttributes attrs{};

  Attribute::write(fp, "format", attrs.format);
  Attribute::write(fp, "format-version", std::int64_t{attrs.format_version});
  assert(attrs.bin_type == BinTable::Type::fixed);
  Attribute::write(fp, "bin-type", std::string{"fixed"});

  auto res_group = fp.createGroup("/resolutions");

  return {fp, chroms, {}, attrs};
}

File MultiResFile::open(std::uint32_t resolution) const {
  const auto match = std::find(resolutions().begin(), resolutions().end(), resolution);

  if (match == resolutions().end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("file \"{}\" does not contain interactions for resolution {}"),
                    path(), resolution));
  }
  return File(
      RootGroup{(*_root_grp)().getGroup(fmt::format(FMT_STRING("/resolutions/{}"), resolution))});
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

File MultiResFile::copy_resolution(const File& clr) {
  SPDLOG_INFO(FMT_STRING("copying {} resolution from {}"), clr.resolution(), clr.uri());
  auto dest = init_resolution(clr.resolution());

  cooler::utils::copy(clr.uri(), dest);
  _resolutions.push_back(clr.resolution());
  std::sort(_resolutions.begin(), _resolutions.end());
  return open(clr.resolution());
}

RootGroup MultiResFile::init_resolution(std::uint32_t resolution) {
  const auto grp = fmt::format(FMT_STRING("/resolutions/{}"), resolution);
  return RootGroup{(*_root_grp)().createGroup(grp, false)};
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

MultiResFile::operator bool() const noexcept { return !!_root_grp; }

std::string MultiResFile::path() const { return (*_root_grp)().getFile().getName(); }

auto MultiResFile::chromosomes() const noexcept -> const Reference& { return _chroms; }

const std::vector<balancing::Method>& MultiResFile::avail_normalizations(
    std::string_view policy) const {
  if (_normalizations.has_value() && _normalizations->first == policy) {
    return _normalizations->second;
  }

  if (policy != "union" && policy != "intersection") {
    throw std::invalid_argument(R"(policy should be either "union" or "intersection")");
  }

  if (_resolutions.empty()) {
    _normalizations.emplace(std::string{policy}, std::vector<balancing::Method>{});
  } else if (_resolutions.size() == 1) {
    _normalizations.emplace(std::string{policy}, open(_resolutions.front()).avail_normalizations());
    return _normalizations->second;
  } else if (policy == "union") {
    _normalizations.emplace(std::string{policy}, avail_normalizations_union());
  } else {
    assert(policy == "intersection");
    _normalizations.emplace(std::string{policy}, avail_normalizations_intersection());
  }
  return _normalizations->second;
}

std::uint32_t MultiResFile::compute_base_resolution(const std::vector<std::uint32_t>& resolutions,
                                                    std::uint32_t target_res) {
  assert(!resolutions.empty());
  const auto base_resolution = resolutions.front();

  if (base_resolution > target_res || target_res % base_resolution != 0) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("resolution {} is not a multiple of base resolution {}"), target_res,
                    base_resolution));
  }

  return *std::find_if(resolutions.rbegin(), resolutions.rend(),
                       [&](const auto res) { return res <= target_res && target_res % res == 0; });
}

HighFive::File MultiResFile::file_handle() {
  assert(_root_grp);
  return (*_root_grp)().getFile();
}
const HighFive::File& MultiResFile::file_handle() const {
  assert(_root_grp);
  return (*_root_grp)().getFile();
}

std::vector<std::uint32_t> MultiResFile::read_resolutions(const HighFive::File& f) {
  [[maybe_unused]] const HighFive::SilenceHDF5 silencer{};  // NOLINT
  try {
    auto root_grp = f.getGroup("/resolutions");

    const auto resolutions_ = root_grp.listObjectNames();
    std::vector<std::uint32_t> resolutions(resolutions_.size());
    std::transform(resolutions_.begin(), resolutions_.end(), resolutions.begin(),
                   [](const auto& res) {
                     return hictk::internal::parse_numeric_or_throw<std::uint32_t>(res);
                   });

    std::sort(resolutions.begin(), resolutions.end());
    return resolutions;
  } catch (const std::exception& e) {
    throw std::runtime_error(fmt::format(FMT_STRING("failed to read resolutions from \"{}\": {}"),
                                         f.getPath(), e.what()));
  }
}

MultiResAttributes MultiResFile::read_attributes(const HighFive::File& f) {
  auto read_or_throw = [&](const auto& key, auto& buff) {
    using T = remove_cvref_t<decltype(buff)>;
    try {
      buff = Attribute::read<T>(f, key);
    } catch (const std::exception& e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("Failed to read attribute \"{}\" from path \"{}\". Reason: {}"),
                      key, f.getPath(), e.what()));
    }
  };

  // Read mandatory attributes
  MultiResAttributes attrs{};
  read_or_throw("format-version", attrs.format_version);
  read_or_throw("format", attrs.format);

  attrs.bin_type = BinTable::Type::fixed;
  if (f.hasAttribute("bin-type")) {
    attrs.bin_type = Attribute::read<std::string>(f, "bin-type") == "fixed"
                         ? BinTable::Type::fixed
                         : BinTable::Type::variable;
  }

  return attrs;
}

std::vector<balancing::Method> MultiResFile::avail_normalizations_union() const {
  assert(_resolutions.size() > 1);

  phmap::flat_hash_set<balancing::Method> norms;
  for (const auto& res : _resolutions) {
    for (const auto& norm : open(res).avail_normalizations()) {
      norms.emplace(norm);
    }
  }

  std::vector norms_sorted(std::make_move_iterator(norms.begin()),
                           std::make_move_iterator(norms.end()));
  std::sort(norms_sorted.begin(), norms_sorted.end());
  return norms_sorted;
}

std::vector<balancing::Method> MultiResFile::avail_normalizations_intersection() const {
  assert(_resolutions.size() > 1);

  phmap::flat_hash_map<balancing::Method, std::uint32_t> norms;
  for (const auto& res : _resolutions) {
    for (const auto& norm : open(res).avail_normalizations()) {
      auto [it, inserted] = norms.try_emplace(norm, std::uint32_t{1});
      if (!inserted) {
        it->second++;
      }
    }
  }

  std::vector<balancing::Method> filtered_norms{};
  for (const auto& [norm, count] : norms) {
    if (count == _resolutions.size()) {
      filtered_norms.emplace_back(norm);
    }
  }

  std::sort(filtered_norms.begin(), filtered_norms.end());
  return filtered_norms;
}

}  // namespace hictk::cooler
