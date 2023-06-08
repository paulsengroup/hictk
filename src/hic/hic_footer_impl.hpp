// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <utility>

#include "hictk/hic/common.hpp"

namespace hictk::internal {

constexpr HiCFooterMetadata::operator bool() const noexcept { return fileOffset >= 0; }

inline bool HiCFooterMetadata::operator==(const HiCFooterMetadata &other) const noexcept {
  return url == other.url && matrix_type == other.matrix_type &&
         normalization == other.normalization && unit == other.unit &&
         resolution == other.resolution && chrom1 == other.chrom1 && chrom2 == other.chrom2;
}

inline bool HiCFooterMetadata::operator!=(const HiCFooterMetadata &other) const noexcept {
  return !(*this == other);
}

inline HiCFooter::HiCFooter(HiCFooterMetadata metadata_) noexcept
    : _metadata(std::move(metadata_)) {}

constexpr HiCFooter::operator bool() const noexcept { return !metadata(); }
inline bool HiCFooter::operator==(const HiCFooter &other) const noexcept {
  return metadata() == other.metadata();
}
inline bool HiCFooter::operator!=(const HiCFooter &other) const noexcept {
  return !(*this == other);
}
constexpr const HiCFooterMetadata &HiCFooter::metadata() const noexcept { return _metadata; }
constexpr HiCFooterMetadata &HiCFooter::metadata() noexcept { return _metadata; }
constexpr const std::string &HiCFooter::url() const noexcept { return metadata().url; }
constexpr MatrixType HiCFooter::matrix_type() const noexcept { return metadata().matrix_type; }
constexpr NormalizationMethod HiCFooter::normalization() const noexcept {
  return metadata().normalization;
}
constexpr MatrixUnit HiCFooter::unit() const noexcept { return metadata().unit; }
constexpr std::int64_t HiCFooter::resolution() const noexcept {
  return std::int64_t(metadata().resolution);
}
constexpr const chromosome &HiCFooter::chrom1() const noexcept { return metadata().chrom1; }
constexpr const chromosome &HiCFooter::chrom2() const noexcept { return metadata().chrom2; }
constexpr std::int64_t HiCFooter::fileOffset() const noexcept { return metadata().fileOffset; }

constexpr const std::vector<double> &HiCFooter::expectedValues() const noexcept {
  return _expectedValues;
}

constexpr const std::vector<double> &HiCFooter::c1Norm() const noexcept { return _c1Norm; }

constexpr const std::vector<double> &HiCFooter::c2Norm() const noexcept {
  if (chrom1() == chrom2()) {
    return _c1Norm;
  }
  return _c2Norm;
}

constexpr std::vector<double> &HiCFooter::expectedValues() noexcept { return _expectedValues; }

constexpr std::vector<double> &HiCFooter::c1Norm() noexcept { return _c1Norm; }

constexpr std::vector<double> &HiCFooter::c2Norm() noexcept {
  if (chrom1() == chrom2()) {
    return _c1Norm;
  }
  return _c2Norm;
}

}  // namespace hictk::internal

template <>
struct std::hash<hictk::internal::HiCFooterMetadata> {
  inline std::size_t operator()(hictk::internal::HiCFooterMetadata const &m) const noexcept {
    return hictk::internal::hash_combine(0, m.url, m.matrix_type, m.normalization, m.unit,
                                         m.resolution, m.chrom1, m.chrom2);
  }
};
