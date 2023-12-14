// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

namespace hictk::cooler {

// Magic values
inline constexpr std::string_view COOL_MAGIC{"HDF5::Cooler"};
inline constexpr std::string_view MCOOL_MAGIC{"HDF5::MCOOL"};
inline constexpr std::string_view SCOOL_MAGIC{"HDF5::SCOOL"};

// clang-format off
inline constexpr std::array<std::string_view, 4> MANDATORY_GROUP_NAMES{
"chroms",
"bins",
"pixels",
"indexes"
};

inline constexpr std::array<std::string_view, 10> MANDATORY_DATASET_NAMES{
"chroms/name",
"chroms/length",
"bins/chrom",
"bins/start",
"bins/end",
"pixels/bin1_id",
"pixels/bin2_id",
"pixels/count",
"indexes/bin1_offset",
"indexes/chrom_offset"
};
// clang-format on

inline constexpr std::uint_fast8_t DEFAULT_COMPRESSION_LEVEL = 6;
inline constexpr std::size_t DEFAULT_HDF5_CHUNK_SIZE = 64ULL << 10U;  // 64KB
inline constexpr double DEFAULT_HDF5_CACHE_W0 = 0.75;
inline constexpr std::size_t DEFAULT_HDF5_DATASET_CACHE_SIZE = 1ULL << 20U;        // 1MB
inline constexpr std::size_t DEFAULT_HDF5_PIXEL_DATASET_CACHE_SIZE = 4ULL << 20U;  // 4MB
inline constexpr std::size_t DEFAULT_HDF5_CACHE_SIZE =                             // 19MB
    (3 * DEFAULT_HDF5_PIXEL_DATASET_CACHE_SIZE) +
    ((MANDATORY_DATASET_NAMES.size() - 3) * DEFAULT_HDF5_DATASET_CACHE_SIZE);

namespace internal {
inline constexpr std::string_view SENTINEL_ATTR_NAME{"format-version"};
inline constexpr std::uint8_t SENTINEL_ATTR_VALUE{255};
}  // namespace internal

}  // namespace hictk::cooler
