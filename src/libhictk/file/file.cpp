// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/file.hpp"

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/pixel_selector.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk {

// NOLINTNEXTLINE(bugprone-exception-escape)
const PixelCoordinates& PixelSelector::coord1() const noexcept {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) -> const PixelCoordinates& {
        using T = std::decay_t<decltype(sel)>;
        if constexpr (std::is_same_v<hic::PixelSelectorAll, T>) {
          static const PixelCoordinates coords{};
          return coords;
        } else {
          return sel.coord1();
        }
      },
      _sel);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
const PixelCoordinates& PixelSelector::coord2() const noexcept {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) -> const PixelCoordinates& {
        using T = std::decay_t<decltype(sel)>;
        if constexpr (std::is_same_v<hic::PixelSelectorAll, T>) {
          static const PixelCoordinates coords{};
          return coords;
        } else {
          return sel.coord2();
        }
      },
      _sel);
}

std::uint64_t PixelSelector::size(bool upper_triangle) const {
  assert(!_sel.valueless_by_exception());
  return std::visit([&](const auto& sel) { return sel.size(upper_triangle); }, _sel);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
const BinTable& PixelSelector::bins() const noexcept {
  assert(!_sel.valueless_by_exception());
  return std::visit([&](const auto& sel) -> const BinTable& { return sel.bins(); }, _sel);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
std::shared_ptr<const BinTable> PixelSelector::bins_ptr() const noexcept {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) -> std::shared_ptr<const BinTable> { return sel.bins_ptr(); }, _sel);
}

PixelSelector PixelSelector::fetch(PixelCoordinates coord1_, PixelCoordinates coord2_) const {
  assert(!_sel.valueless_by_exception());
  return std::visit(
      [&](const auto& sel) -> PixelSelector {
        using T = remove_cvref_t<decltype(sel)>;
        if constexpr (std::is_same_v<hic::PixelSelectorAll, T>) {
          throw std::runtime_error(
              "calling fetch() on a PixelSelector instance set up to fetch genome-wide matrices "
              "from .hic files is not supported");
        } else {
          return PixelSelector{sel.fetch(coord1_, coord2_), _weights};
        }
      },
      _sel);
}

const balancing::Weights& PixelSelector::weights() const noexcept {
  assert(_weights);
  return *_weights;
}

[[nodiscard]] static cooler::File open_single_res_cooler(std::string_view uri,
                                                         std::optional<std::uint32_t> resolution_) {
  cooler::File f{uri};
  if (resolution_.value_or(f.resolution()) != f.resolution()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("found an unexpected resolution while opening file at \"{}\": "
                               "expected {}, found {}."),
                    uri, *resolution_, f.resolution()));  // NOLINT(*-unchecked-optional-access)
  }
  return f;
}

[[nodiscard]] static std::uint32_t try_infer_resolution(std::string_view uri) {
  const auto resolutions = [&]() {
    try {
      return cooler::utils::list_resolutions(uri, false);
    } catch (const std::exception&) {
      return std::vector<std::uint32_t>{};
    }
  }();

  if (resolutions.size() == 1) {
    return resolutions.front();
  }
  throw std::runtime_error(
      "resolution is required when opening .mcool files with more than one resolution.");
}

[[noreturn]] static void raise_invalid_resolution_except(
    const cooler::utils::ValidationStatusCooler& status, std::string_view path,
    std::string_view uri, std::uint32_t resolution) {
  const auto grp = fmt::format(FMT_STRING("resolutions/{}"), resolution);
  const auto resolution_is_missing = status.missing_groups.size() == 1 &&
                                     status.missing_groups.front().find(grp) != std::string::npos;
  if (resolution_is_missing) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find resolution {} in file \"{}\""), resolution, path));
  }
  throw std::runtime_error(fmt::format(FMT_STRING("\"{}\" does not look like a valid Cooler file:\n"
                                                  "Validation report:\n{}"),
                                       uri, status));
}

[[nodiscard]] static std::variant<cooler::File, hic::File> file_ctor_helper(
    std::string_view uri, std::optional<std::uint32_t> resolution_, hic::MatrixType type,
    hic::MatrixUnit unit) {
  const auto [path, grp] = cooler::parse_cooler_uri(uri);
  if (hic::utils::is_hic_file(path)) {
    return {hic::File(path, resolution_, type, unit)};
  }

  if (type != hic::MatrixType::observed) {
    throw std::runtime_error(
        "matrix type should always be \"observed\" when reading Cooler files.");
  }

  if (unit != hic::MatrixUnit::BP) {
    throw std::runtime_error("matrix unit should always be \"BP\" when reading Cooler files.");
  }

  if (cooler::utils::is_cooler(uri)) {
    return {open_single_res_cooler(uri, resolution_)};
  }

  if (!resolution_.has_value()) {
    resolution_ = try_infer_resolution(uri);
  }

  const auto new_uri = fmt::format(FMT_STRING("{}::/resolutions/{}"), path, *resolution_);

  const auto status = cooler::utils::is_cooler(new_uri);
  if (status) {
    return {cooler::File(new_uri, cooler::DEFAULT_HDF5_CACHE_SIZE * 4, false)};
  }

  raise_invalid_resolution_except(status, path, new_uri, *resolution_);
}

File::File(cooler::File clr) : _fp(std::move(clr)) {}
File::File(hic::File hf) : _fp(std::move(hf)) {}
File::File(std::string_view uri, std::optional<std::uint32_t> resolution_, hic::MatrixType type,
           hic::MatrixUnit unit)
    : _fp(file_ctor_helper(uri, resolution_, type, unit)) {}

std::string File::uri() const {
  assert(!_fp.valueless_by_exception());
  // NOLINTBEGIN(bugprone-branch-clone)
  return std::visit(
      [&](const auto& fp) {
        using T = std::decay_t<decltype(fp)>;
        if constexpr (std::is_same_v<hic::File, T>) {
          return fp.path();
        } else {
          return fp.uri();
        }
      },
      _fp);
  // NOLINTEND(bugprone-branch-clone)
}

std::string File::path() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) { return fp.path(); }, _fp);
}

auto File::chromosomes() const -> const Reference& {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) -> const Reference& { return fp.chromosomes(); }, _fp);
}

auto File::bins() const -> const BinTable& {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) -> const BinTable& { return fp.bins(); }, _fp);
}

std::shared_ptr<const BinTable> File::bins_ptr() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& f) -> std::shared_ptr<const BinTable> { return f.bins_ptr(); },
                    _fp);
}

std::uint32_t File::resolution() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) { return fp.resolution(); }, _fp);
}

std::uint64_t File::nbins() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) { return fp.nbins(); }, _fp);
}

std::uint64_t File::nchroms(bool include_ALL) const {
  assert(!_fp.valueless_by_exception());

  return std::visit(
      [&](const auto& fp) {
        using T = remove_cvref_t<decltype(fp)>;
        if constexpr (std::is_same_v<hic::File, T>) {
          return fp.nchroms(include_ALL);
        } else {
          return fp.nchroms();
        }
      },
      _fp);
}

PixelSelector File::fetch(const balancing::Method& normalization) const {
  assert(!_fp.valueless_by_exception());
  return std::visit(
      [&](const auto& fp) {
        return PixelSelector{fp.fetch(normalization), fp.normalization_ptr(normalization)};
      },
      _fp);
}

PixelSelector File::fetch(std::string_view range, const balancing::Method& normalization,
                          QUERY_TYPE query_type) const {
  assert(!_fp.valueless_by_exception());
  return fetch(range, range, normalization, query_type);
}

PixelSelector File::fetch(std::string_view chrom_name, std::uint32_t start, std::uint32_t end,
                          const balancing::Method& normalization) const {
  return fetch(chrom_name, start, end, chrom_name, start, end, normalization);
}

PixelSelector File::fetch(std::string_view range1, std::string_view range2,
                          const balancing::Method& normalization, QUERY_TYPE query_type) const {
  assert(!_fp.valueless_by_exception());
  return std::visit(
      [&](const auto& fp) {
        return PixelSelector{fp.fetch(range1, range2, normalization, query_type),
                             fp.normalization_ptr(normalization)};
      },
      _fp);
}

PixelSelector File::fetch(std::string_view chrom1_name, std::uint32_t start1, std::uint32_t end1,
                          std::string_view chrom2_name, std::uint32_t start2, std::uint32_t end2,
                          const balancing::Method& normalization) const {
  assert(!_fp.valueless_by_exception());
  return std::visit(
      [&](const auto& fp) {
        return PixelSelector{
            fp.fetch(chrom1_name, start1, end1, chrom2_name, start2, end2, normalization),
            fp.normalization_ptr(normalization)};
      },
      _fp);
}

bool File::has_normalization(std::string_view normalization) const {
  assert(!_fp.valueless_by_exception());
  return std::visit([&](const auto& fp) { return fp.has_normalization(normalization); }, _fp);
}
std::vector<balancing::Method> File::avail_normalizations() const {
  assert(!_fp.valueless_by_exception());
  return std::visit([](const auto& fp) { return fp.avail_normalizations(); }, _fp);
}

const balancing::Weights& File::normalization(std::string_view normalization_) const {
  assert(!_fp.valueless_by_exception());
  if (std::holds_alternative<cooler::File>(_fp)) {
    return std::get<cooler::File>(_fp).normalization(normalization_);
  }
  return std::get<hic::File>(_fp).normalization(normalization_);
}

}  // namespace hictk
