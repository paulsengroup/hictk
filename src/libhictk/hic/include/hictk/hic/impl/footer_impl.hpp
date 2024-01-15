// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/hash.hpp"
#include "hictk/hic/common.hpp"

namespace hictk::hic::internal {

constexpr HiCFooterMetadata::operator bool() const noexcept { return matrixMetadataOffset >= 0; }

inline bool HiCFooterMetadata::operator==(const HiCFooterMetadata &other) const noexcept {
  return url == other.url && matrix_type == other.matrix_type &&
         normalization == other.normalization && unit == other.unit &&
         resolution == other.resolution && chrom1 == other.chrom1 && chrom2 == other.chrom2;
}

inline bool HiCFooterMetadata::operator!=(const HiCFooterMetadata &other) const noexcept {
  return !(*this == other);
}

inline HiCFooter::HiCFooter(Index index_, HiCFooterMetadata metadata_,
                            std::vector<double> expected_values,
                            std::shared_ptr<balancing::Weights> weights1,
                            std::shared_ptr<balancing::Weights> weights2) noexcept
    : _index(std::move(index_)),
      _metadata(std::move(metadata_)),
      _expectedValues(std::move(expected_values)),
      _weights1(std::move(weights1)),
      _weights2(std::move(weights2)) {}

constexpr HiCFooter::operator bool() const noexcept { return !metadata(); }
inline bool HiCFooter::operator==(const HiCFooter &other) const noexcept {
  return metadata() == other.metadata();
}
inline bool HiCFooter::operator!=(const HiCFooter &other) const noexcept {
  return !(*this == other);
}
constexpr const HiCFooterMetadata &HiCFooter::metadata() const noexcept { return _metadata; }
constexpr HiCFooterMetadata &HiCFooter::metadata() noexcept { return _metadata; }
inline const Index &HiCFooter::index() const noexcept { return _index; }
constexpr const std::string &HiCFooter::url() const noexcept { return metadata().url; }
constexpr MatrixType HiCFooter::matrix_type() const noexcept { return metadata().matrix_type; }
inline balancing::Method HiCFooter::normalization() const noexcept {
  return metadata().normalization;
}
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

inline const balancing::Weights &HiCFooter::weights1() const noexcept {
  assert(_weights1);
  return *_weights1;
}

inline const balancing::Weights &HiCFooter::weights2() const noexcept {
  assert(_weights2);
  return *_weights2;
}

constexpr std::vector<double> &HiCFooter::expectedValues() noexcept { return _expectedValues; }
}  // namespace hictk::hic::internal

inline std::size_t std::hash<hictk::hic::internal::HiCFooterMetadata>::operator()(
    hictk::hic::internal::HiCFooterMetadata const &m) const noexcept {
  return hictk::internal::hash_combine(0, m.url, m.matrix_type, m.normalization.to_string(), m.unit,
                                       m.resolution, m.chrom1, m.chrom2);
}

inline std::size_t std::hash<hictk::hic::internal::HiCFooter>::operator()(
    hictk::hic::internal::HiCFooter const &f) const noexcept {
  return std::hash<hictk::hic::internal::HiCFooterMetadata>{}(f.metadata());
}
