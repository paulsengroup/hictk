// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"

namespace hictk::hic::internal {

constexpr HiCFooterMetadata::operator bool() const noexcept { return matrixMetadataOffset >= 0; }

constexpr HiCFooter::operator bool() const noexcept { return !metadata(); }

constexpr const HiCFooterMetadata &HiCFooter::metadata() const noexcept { return _metadata; }

constexpr HiCFooterMetadata &HiCFooter::metadata() noexcept { return _metadata; }

constexpr const std::string &HiCFooter::path() const noexcept { return metadata().url; }

constexpr MatrixType HiCFooter::matrix_type() const noexcept { return metadata().matrix_type; }

constexpr MatrixUnit HiCFooter::unit() const noexcept { return metadata().unit; }

constexpr std::uint32_t HiCFooter::resolution() const noexcept { return metadata().resolution; }

constexpr const Chromosome &HiCFooter::chrom1() const noexcept { return metadata().chrom1; }

constexpr const Chromosome &HiCFooter::chrom2() const noexcept { return metadata().chrom2; }

constexpr std::int64_t HiCFooter::fileOffset() const noexcept {
  return metadata().matrixMetadataOffset;
}

constexpr const std::vector<double> &HiCFooter::expectedValues() const noexcept {
  return _expectedValues;
}

constexpr std::vector<double> &HiCFooter::expectedValues() noexcept { return _expectedValues; }
}  // namespace hictk::hic::internal
