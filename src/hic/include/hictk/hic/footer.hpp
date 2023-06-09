// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/index.hpp"

namespace hictk::hic::internal {
struct HiCFooterMetadata {
  std::string url{};
  MatrixType matrix_type{MatrixType::observed};
  NormalizationMethod normalization{NormalizationMethod::NONE};
  MatrixUnit unit{MatrixUnit::BP};
  std::uint32_t resolution{std::numeric_limits<std::uint32_t>::max()};
  Chromosome chrom1{};
  Chromosome chrom2{};
  std::int64_t fileOffset{-1};

  constexpr explicit operator bool() const noexcept;
  bool operator==(const HiCFooterMetadata &other) const noexcept;
  bool operator!=(const HiCFooterMetadata &other) const noexcept;
};

class HiCFooter {
  Index _index{};
  HiCFooterMetadata _metadata{};
  std::vector<double> _expectedValues{};
  std::vector<double> _c1Norm{};
  std::vector<double> _c2Norm{};

 public:
  HiCFooter() = default;
  explicit HiCFooter(Index index_, HiCFooterMetadata metadata_) noexcept;

  constexpr explicit operator bool() const noexcept;
  bool operator==(const HiCFooter &other) const noexcept;
  bool operator!=(const HiCFooter &other) const noexcept;

  [[nodiscard]] constexpr const HiCFooterMetadata &metadata() const noexcept;
  [[nodiscard]] constexpr HiCFooterMetadata &metadata() noexcept;
  [[nodiscard]] const Index &index() const noexcept;

  [[nodiscard]] constexpr const std::string &url() const noexcept;
  [[nodiscard]] constexpr MatrixType matrix_type() const noexcept;
  [[nodiscard]] constexpr NormalizationMethod normalization() const noexcept;
  [[nodiscard]] constexpr MatrixUnit unit() const noexcept;
  [[nodiscard]] constexpr std::uint32_t resolution() const noexcept;
  [[nodiscard]] constexpr const Chromosome &chrom1() const noexcept;
  [[nodiscard]] constexpr const Chromosome &chrom2() const noexcept;
  [[nodiscard]] constexpr std::int64_t fileOffset() const noexcept;

  [[nodiscard]] constexpr const std::vector<double> &expectedValues() const noexcept;
  [[nodiscard]] constexpr const std::vector<double> &c1Norm() const noexcept;
  [[nodiscard]] constexpr const std::vector<double> &c2Norm() const noexcept;

  [[nodiscard]] constexpr std::vector<double> &expectedValues() noexcept;
  [[nodiscard]] constexpr std::vector<double> &c1Norm() noexcept;
  [[nodiscard]] constexpr std::vector<double> &c2Norm() noexcept;
};
}  // namespace hictk::hic::internal

template <>
struct std::hash<hictk::hic::internal::HiCFooterMetadata> {
  inline std::size_t operator()(hictk::hic::internal::HiCFooterMetadata const &m) const noexcept;
};

template <>
struct std::hash<hictk::hic::internal::HiCFooter> {
  inline std::size_t operator()(hictk::hic::internal::HiCFooter const &f) const noexcept;
};
#include "../../../footer_impl.hpp"
