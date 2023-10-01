// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hictk/cooler/attribute.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/transformers/coarsen.hpp"

namespace hictk::cooler {

inline bool MultiResAttributes::operator==(const MultiResAttributes& other) const noexcept {
  return format == other.format && format_version == other.format_version &&
         bin_type == other.bin_type;
}
inline bool MultiResAttributes::operator!=(const MultiResAttributes& other) const noexcept {
  return !(*this == other);
}

inline MultiResFile::MultiResFile(HighFive::File fp, Reference chroms,
                                  std::vector<std::uint32_t> resolutions, MultiResAttributes attrs)
    : _root_grp(std::make_unique<RootGroup>(RootGroup{fp.getGroup("/")})),
      _resolutions(std::move(resolutions)),
      _attrs(std::move(attrs)),
      _chroms(std::move(chroms)) {
  if (_chroms.empty() && !_resolutions.empty()) {
    _chroms = open(_resolutions.back()).chromosomes();
  }
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

inline MultiResFile::MultiResFile(const std::filesystem::path& path, unsigned int mode)
    : MultiResFile(HighFive::File(path.string(), mode), {},
                   read_resolutions(HighFive::File(path.string())),
                   read_attributes(HighFive::File(path.string()))) {}

inline MultiResFile MultiResFile::create(const std::filesystem::path& path, const Reference& chroms,
                                         bool force_overwrite) {
  if (!force_overwrite && std::filesystem::exists(path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to initialize file \"{}\": file already exists"), path));
  }

  if (force_overwrite) {
    std::filesystem::remove(path);
  }
  HighFive::File fp(path.string(), HighFive::File::Create);

  MultiResAttributes attrs{};

  Attribute::write(fp, "format", attrs.format);
  Attribute::write(fp, "format-version", std::int64_t(attrs.format_version));
  if (attrs.bin_type) {
    Attribute::write(fp, "bin-type", *attrs.bin_type);
  }

  auto res_group = fp.createGroup("/resolutions");

  return {fp, chroms, {}, attrs};
}

template <typename ResolutionIt>
inline MultiResFile MultiResFile::create(const std::filesystem::path& path, const File& base,
                                         ResolutionIt first_res, ResolutionIt last_res,
                                         bool force_overwrite) {
  std::vector<std::uint32_t> resolutions_{first_res, last_res};
  std::sort(resolutions_.begin(), resolutions_.end());
  const auto base_res = base.bin_size();
  const auto tgt_res = resolutions_.front();

  for (const auto& res : resolutions_) {
    if (base_res > tgt_res || res % base_res != 0) {
      throw std::runtime_error(fmt::format(
          FMT_STRING("resolution {} is not a multiple of base resolution {}"), res, base_res));
    }
  }

  auto mclr = MultiResFile::create(path, base.chromosomes(), force_overwrite);
  SPDLOG_INFO(FMT_STRING("Copying {} resolution from \"{}\""), base_res, base.path());
  mclr.copy_resolution(base);

  std::vector<ThinPixel<std::int32_t>> buffer{500'000};
  for (std::size_t i = 1; i < resolutions_.size(); ++i) {
    const auto tgt_resolution = resolutions_[i];
    const auto base_resolution = compute_base_resolution(mclr._resolutions, tgt_resolution);

    auto attributes = base.has_float_pixels() ? Attributes::init<double>(tgt_resolution)
                                              : Attributes::init<std::int32_t>(tgt_resolution);
    attributes.assembly = base.attributes().assembly;
    auto clr =
        base.has_float_pixels()
            ? File::create<double>(mclr.init_resolution(tgt_resolution), base.chromosomes(),
                                   tgt_resolution, attributes, DEFAULT_HDF5_CACHE_SIZE)
            : File::create<std::int32_t>(mclr.init_resolution(tgt_resolution), base.chromosomes(),
                                         tgt_resolution, attributes, DEFAULT_HDF5_CACHE_SIZE);
    SPDLOG_INFO(FMT_STRING("Generating {} resolution from {} ({}x)"), tgt_resolution,
                base_resolution, tgt_resolution / base_resolution);

    MultiResFile::coarsen(mclr.open(base_res), clr, buffer);
    mclr._resolutions.push_back(tgt_resolution);
  }

  return mclr;
}

constexpr const std::vector<std::uint32_t>& MultiResFile::resolutions() const noexcept {
  return _resolutions;
}

constexpr const MultiResAttributes& MultiResFile::attributes() const noexcept { return _attrs; }

inline File MultiResFile::open(std::uint32_t resolution) const {
  return File(
      RootGroup{(*_root_grp)().getGroup(fmt::format(FMT_STRING("/resolutions/{}"), resolution))});
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

inline File MultiResFile::copy_resolution(const File& clr) {
  SPDLOG_INFO(FMT_STRING("copying {} resolution from {}"), clr.bin_size(), clr.uri());
  auto dest = init_resolution(clr.bin_size());

  cooler::utils::copy(clr.uri(), dest);
  _resolutions.push_back(clr.bin_size());
  std::sort(_resolutions.begin(), _resolutions.end());
  return open(clr.bin_size());
}

template <typename N>
inline File MultiResFile::create_resolution(std::uint32_t resolution, Attributes attributes) {
  const auto base_resolution = compute_base_resolution(resolutions(), resolution);

  std::vector<ThinPixel<std::int32_t>> buffer{500'000};
  auto base_clr = open(base_resolution);
  attributes.assembly = base_clr.attributes().assembly;
  attributes.bin_size = resolution;
  {
    auto clr = File::create<N>(init_resolution(resolution), base_clr.chromosomes(), resolution,
                               attributes, DEFAULT_HDF5_CACHE_SIZE);
    MultiResFile::coarsen(base_clr, clr, buffer);
  }

  _resolutions.push_back(resolution);
  std::sort(_resolutions.begin(), _resolutions.end());

  return open(resolution);
}

inline RootGroup MultiResFile::init_resolution(std::uint32_t resolution) {
  const auto grp = fmt::format(FMT_STRING("/resolutions/{}"), resolution);
  return RootGroup{(*_root_grp)().createGroup(grp, false)};
}  // NOLINT(clang-analyzer-cplusplus.NewDeleteLeaks)

inline MultiResFile::operator bool() const noexcept { return !!_root_grp; }

inline std::string MultiResFile::path() const { return (*_root_grp)().getFile().getName(); }

inline auto MultiResFile::chromosomes() const noexcept -> const Reference& {
  return _chroms;
}

[[nodiscard]] inline std::uint32_t MultiResFile::compute_base_resolution(
    const std::vector<std::uint32_t>& resolutions, std::uint32_t target_res) {
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

inline void MultiResFile::coarsen(const File& clr1, File& clr2,
                                  std::vector<ThinPixel<std::int32_t>>& buffer) {
  SPDLOG_INFO(FMT_STRING("generating {} resolution from {} ({}x)"), clr2.bin_size(),
              clr1.bin_size(), clr2.bin_size() / clr1.bin_size());
  auto sel1 = clr1.fetch();
  auto sel2 = transformers::CoarsenPixels(sel1.begin<std::int32_t>(), sel1.end<std::int32_t>(),
                                          clr1.bins_ptr(), clr2.bin_size() / clr1.bin_size());

  const auto update_frequency =
      std::max(std::size_t(1'000'000), (clr1.dataset("pixels/bin1_id").size() / 100));

  auto first = sel2.begin();
  auto last = sel2.end();
  buffer.clear();

  auto t0 = std::chrono::steady_clock::now();
  for (std::size_t j = 0; first != last; ++j) {
    buffer.emplace_back(*first);
    if (buffer.size() == buffer.capacity()) {
      clr2.append_pixels(buffer.begin(), buffer.end());
      buffer.clear();
    }
    if (j == update_frequency) {
      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      const auto bin1 = clr2.bins().at(first->bin1_id);
      SPDLOG_INFO(FMT_STRING("[{} -> {}] processing {:ucsc} at {:.0f} pixels/s..."),
                  clr1.bin_size(), clr2.bin_size(), bin1, double(update_frequency) / delta);
      t0 = t1;
      j = 0;
    }
    ++first;
  }
  if (!buffer.empty()) {
    clr2.append_pixels(buffer.begin(), buffer.end());
  }
}

inline std::vector<std::uint32_t> MultiResFile::read_resolutions(const HighFive::File& f) {
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

inline MultiResAttributes MultiResFile::read_attributes(const HighFive::File& f) {
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

  if (f.hasAttribute("bin-type")) {
    attrs.bin_type = Attribute::read<std::string>(f, "bin-type");
  } else {
    attrs.bin_type.reset();
  }

  return attrs;
}

}  // namespace hictk::cooler
