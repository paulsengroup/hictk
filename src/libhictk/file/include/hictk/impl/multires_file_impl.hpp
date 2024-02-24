// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/cooler/utils.hpp"
#include "hictk/cooler/validation.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"

namespace hictk {

inline MultiResFile::MultiResFile(const cooler::MultiResFile& mclr)
    : MultiResFile(mclr.path(), hic::MatrixType::observed, hic::MatrixUnit::BP) {}

inline MultiResFile::MultiResFile(const hic::File& hf)
    : MultiResFile(hf.path(), hf.matrix_type(), hf.matrix_unit()) {}

inline MultiResFile::MultiResFile(std::string uri, hic::MatrixType type_, hic::MatrixUnit unit_)
    : _path(std::move(uri)), _type(type_), _unit(unit_) {
  if (hic::utils::is_hic_file(_path)) {
    _resolutions = hic::utils::list_resolutions(_path);
    const hic::File hf(_path, _resolutions.back());
    _chroms = hf.chromosomes();
    _format = "HIC";
    _format_version = static_cast<std::uint8_t>(hf.version());
    _bin_type = "fixed";
    return;
  }
  if (!cooler::utils::is_multires_file(_path)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("file is not in .hic or .mcool format: \"{}\""), _path));
  }

  if (type_ != hic::MatrixType::observed) {
    throw std::runtime_error(
        "matrix type should always be \"observed\" when opening .mcool files.");
  }

  if (unit_ != hic::MatrixUnit::BP) {
    throw std::runtime_error("matrix unit should always be \"BP\" when opening .mcool files.");
  }

  const auto mclr = cooler::MultiResFile(_path);
  _chroms = mclr.chromosomes();
  _resolutions = mclr.resolutions();
  _format = cooler::MCOOL_MAGIC;
  _format_version = mclr.attributes().format_version;
  _bin_type = mclr.attributes().bin_type.has_value() ? *mclr.attributes().bin_type : "fixed";
}

inline std::string MultiResFile::path() const { return _path; }

inline bool MultiResFile::is_hic() const noexcept { return _format == "HIC"; }
inline bool MultiResFile::is_mcool() const noexcept { return !is_hic(); }

constexpr hic::MatrixType MultiResFile::matrix_type() const noexcept { return _type; }
constexpr hic::MatrixUnit MultiResFile::matrix_unit() const noexcept { return _unit; }
inline std::string_view MultiResFile::format() const noexcept { return _format; }
constexpr std::uint8_t MultiResFile::version() const noexcept { return _format_version; }
inline std::string_view MultiResFile::bin_type() const noexcept { return _bin_type; }

constexpr const std::vector<std::uint32_t>& MultiResFile::resolutions() const noexcept {
  return _resolutions;
}

inline const Reference& MultiResFile::chromosomes() const noexcept { return _chroms; }

inline File MultiResFile::open(std::uint32_t resolution) const { return File(_path, resolution); }

}  // namespace hictk
