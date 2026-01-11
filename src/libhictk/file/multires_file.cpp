// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/multires_file.hpp"

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/balancing//methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/utils.hpp"
#include "hictk/hic/validation.hpp"
#include "hictk/reference.hpp"

namespace hictk {

MultiResFile::MultiResFile(const cooler::MultiResFile& mclr)
    : MultiResFile(mclr.path(), hic::MatrixType::observed, hic::MatrixUnit::BP) {}

MultiResFile::MultiResFile(const hic::File& hf)
    : MultiResFile(hf.path(), hf.matrix_type(), hf.matrix_unit()) {}

[[nodiscard]] static MultiResFile::Format infer_format(std::string_view path) {
  if (hic::utils::is_hic_file(path)) {
    return MultiResFile::Format::HIC;
  }

  if (cooler::utils::is_multires_file(path)) {
    return MultiResFile::Format::MCOOL;
  }
  throw std::runtime_error(
      fmt::format(FMT_STRING("file is not in .hic or .mcool format: \"{}\""), path));
}

MultiResFile::MultiResFile(std::string uri, hic::MatrixType type_, hic::MatrixUnit unit_)
    : _path(std::move(uri)), _type(type_), _unit(unit_), _format(infer_format(_path)) {
  if (_format == Format::HIC) {
    _resolutions = hic::utils::list_resolutions(_path);
    const hic::File hf(_path, _resolutions.back());
    _chroms = hf.chromosomes();
    _format_version = static_cast<std::uint8_t>(hf.version());
    return;
  }

  assert(_format == Format::MCOOL);
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
  _format_version = mclr.attributes().format_version;
  _bin_type = mclr.attributes().bin_type;
}

std::string MultiResFile::path() const { return _path; }

const Reference& MultiResFile::chromosomes() const noexcept { return _chroms; }

const std::vector<balancing::Method>& MultiResFile::avail_normalizations(
    std::string_view policy) const {
  if (_normalizations.has_value() && _normalizations->first == policy) {
    return _normalizations->second;
  }

  if (is_hic()) {
    _normalizations.emplace(std::string{policy},
                            hic::utils::list_normalizations(_path, policy, _type, _unit));
  } else {
    assert(is_mcool());
    _normalizations.emplace(std::string{policy},
                            cooler::MultiResFile(_path).avail_normalizations(policy));
  }

  return _normalizations->second;
}

File MultiResFile::open(std::uint32_t resolution) const { return File(_path, resolution); }

}  // namespace hictk
