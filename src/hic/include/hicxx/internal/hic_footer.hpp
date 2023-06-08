// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "hicxx/internal/common.hpp"

namespace hicxx::internal {
struct HiCFooterMetadata {
    std::string url{};
    MatrixType matrix_type{MatrixType::observed};
    NormalizationMethod normalization{NormalizationMethod::NONE};
    MatrixUnit unit{MatrixUnit::BP};
    std::int32_t resolution{-1};
    chromosome chrom1{};
    chromosome chrom2{};
    std::int64_t fileOffset{-1};

    constexpr explicit operator bool() const noexcept;
    bool operator==(const HiCFooterMetadata &other) const noexcept;
    bool operator!=(const HiCFooterMetadata &other) const noexcept;
};

class HiCFooter {
    HiCFooterMetadata _metadata{};
    std::vector<double> _expectedValues{};
    std::vector<double> _c1Norm{};
    std::vector<double> _c2Norm{};

   public:
    HiCFooter() = default;
    explicit HiCFooter(HiCFooterMetadata metadata_) noexcept;

    constexpr explicit operator bool() const noexcept;
    bool operator==(const HiCFooter &other) const noexcept;
    bool operator!=(const HiCFooter &other) const noexcept;

    [[nodiscard]] constexpr const HiCFooterMetadata &metadata() const noexcept;
    [[nodiscard]] constexpr HiCFooterMetadata &metadata() noexcept;

    [[nodiscard]] constexpr const std::string &url() const noexcept;
    [[nodiscard]] constexpr MatrixType matrix_type() const noexcept;
    [[nodiscard]] constexpr NormalizationMethod normalization() const noexcept;
    [[nodiscard]] constexpr MatrixUnit unit() const noexcept;
    [[nodiscard]] constexpr std::int64_t resolution() const noexcept;
    [[nodiscard]] constexpr const chromosome &chrom1() const noexcept;
    [[nodiscard]] constexpr const chromosome &chrom2() const noexcept;
    [[nodiscard]] constexpr std::int64_t fileOffset() const noexcept;

    [[nodiscard]] constexpr const std::vector<double> &expectedValues() const noexcept;
    [[nodiscard]] constexpr const std::vector<double> &c1Norm() const noexcept;
    [[nodiscard]] constexpr const std::vector<double> &c2Norm() const noexcept;

    [[nodiscard]] constexpr std::vector<double> &expectedValues() noexcept;
    [[nodiscard]] constexpr std::vector<double> &c1Norm() noexcept;
    [[nodiscard]] constexpr std::vector<double> &c2Norm() noexcept;
};
}  // namespace hicxx::internal

#include "../../../hic_footer_impl.hpp"
