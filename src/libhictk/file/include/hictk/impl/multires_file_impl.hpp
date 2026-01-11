// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string_view>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/cooler/common.hpp"
#include "hictk/hic/common.hpp"

namespace hictk {

constexpr hic::MatrixType MultiResFile::matrix_type() const noexcept { return _type; }
constexpr hic::MatrixUnit MultiResFile::matrix_unit() const noexcept { return _unit; }
constexpr std::uint8_t MultiResFile::version() const noexcept { return _format_version; }
constexpr BinTable::Type MultiResFile::bin_type() const noexcept { return _bin_type; }

constexpr const std::vector<std::uint32_t>& MultiResFile::resolutions() const noexcept {
  return _resolutions;
}

constexpr bool MultiResFile::is_hic() const noexcept { return _format == Format::HIC; }
constexpr bool MultiResFile::is_mcool() const noexcept { return _format == Format::MCOOL; }

constexpr std::string_view MultiResFile::format() const noexcept {
  switch (_format) {
    case Format::HIC:
      return "HIC";
    case Format::MCOOL:
      return cooler::MCOOL_MAGIC;
    default:
      return "unknown";
  }
}

}  // namespace hictk
