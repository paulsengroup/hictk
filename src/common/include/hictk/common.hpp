// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <stdexcept>
#include <string_view>
#include <utility>

#include "hictk/version.hpp"

namespace hictk {

inline const std::string_view HICTK_VERSION_STRING{hictk::config::version::str()};

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

inline constexpr std::size_t DEFAULT_HDF5_DATASET_ITERATOR_BUFFER_SIZE = 32ULL << 10U;  // 32K

namespace internal {
inline constexpr std::string_view SENTINEL_ATTR_NAME{"format-version"};
inline constexpr std::uint8_t SENTINEL_ATTR_VALUE{255};
}  // namespace internal

[[nodiscard]] constexpr bool ndebug_defined() noexcept {
#ifdef NDEBUG
  return true;
#else
  return false;
#endif
}

[[nodiscard]] constexpr bool ndebug_not_defined() noexcept { return !ndebug_defined(); }

#if defined(__GNUC__) || defined(__builtin_unreachable)
#define HICTK_UNREACHABLE_CODE __builtin_unreachable()
#elif defined(_MSC_VER)
#define HICTK_UNREACHABLE_CODE __assume(0)
#else
#define HICTK_UNREACHABLE_CODE
#endif

[[nodiscard]] constexpr bool noexcept_move_ctor() noexcept {
#if defined(__GNUC__) && !defined(__clang__)
  return __GNUC__ > 7;
#else
  return true;
#endif
}

[[nodiscard]] constexpr bool noexcept_move_assigment_op() noexcept {
#if defined(__GNUC__) && defined(__clang__)
  return __clang_major__ > 8;
#else
  return noexcept_move_ctor();
#endif
}

[[noreturn]] inline void unreachable_code() {
  if constexpr (ndebug_not_defined()) {
    throw std::logic_error("Unreachable code");
  }
  HICTK_UNREACHABLE_CODE;
}

struct identity {
  template <typename T>
  [[nodiscard]] constexpr T &&operator()(T &&a) const noexcept {
    return std::forward<T>(a);
  }
  using is_transparent = void;
};

// to avoid useless casts (see https://github.com/nlohmann/json/issues/2893#issuecomment-889152324)
template <typename T, typename U>
[[maybe_unused]] [[nodiscard]] constexpr T conditional_static_cast(U value) {
  if constexpr (std::is_same_v<T, U>) {
    return value;
  } else {
    return static_cast<T>(value);
  }
}

// helper function to construct unique/shared ptr with a custom deleter fx
template <auto fn>
struct DeleterFromFn {
  template <typename T>
  constexpr void operator()(T *arg) const {
    fn(arg);
  }
};

}  // namespace hictk
