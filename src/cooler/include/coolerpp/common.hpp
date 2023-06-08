// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <array>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>

#include "coolerpp/internal/version.hpp"

namespace coolerpp {

inline const std::string_view COOLERPP_VERSION_STRING{coolerpp::config::version::str()};

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
#define COOLERPP_UNREACHABLE_CODE __builtin_unreachable()
#elif defined(_MSC_VER)
#define COOLERPP_UNREACHABLE_CODE __assume(0)
#else
#define COOLERPP_UNREACHABLE_CODE
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
  COOLERPP_UNREACHABLE_CODE;
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

// metaprogramming stuff

template <typename T>
struct remove_cvref {
  using type = typename std::remove_cv<typename std::remove_reference<T>::type>::type;
};

template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

template <typename T>
struct is_string
    : public std::disjunction<std::is_same<char *, typename std::decay_t<T>>,
                              std::is_same<const char *, typename std::decay_t<T>>,
                              std::is_same<std::string, typename std::decay_t<T>>,
                              std::is_same<std::string_view, typename std::decay_t<T>>> {};

template <typename T>
constexpr bool is_string_v = is_string<T>::value;

template <typename Operation, typename Operand>
struct is_unary_operation : public std::is_invocable<Operation, Operand> {};

template <typename Operation, typename Operand>
constexpr bool is_unary_operation_v = is_unary_operation<Operation, Operand>::value;

namespace internal {
// Adapted from https://stackoverflow.com/a/29634934
// clang-format off
template <typename T>
auto is_iterable_impl(int)
    -> decltype(std::begin(std::declval<T &>()) != std::end(std::declval<T &>()),  // begin/end and operator !=
                void(),                                                            // Handle evil operator ,
                ++std::declval<decltype(std::begin(std::declval<T &>())) &>(),     // operator ++
                void(*std::begin(std::declval<T &>())),                            // operator*
                std::true_type{});
// clang-format on

template <typename T>
std::false_type is_iterable_impl(...);
}  // namespace internal

template <typename T, typename = void>
struct is_iterable : std::false_type {};

// clang-format off
template <typename T>
struct is_iterable<T, std::void_t<decltype(internal::is_iterable_impl<T>(0))>>
    : std::true_type {};
// clang-format on

template <typename T>
constexpr bool is_iterable_v = is_iterable<T>::value;

template <typename Operation, typename It, typename = void>
constexpr bool is_unary_operation_on_iterator = false;

template <typename Operation, typename It>
constexpr bool is_unary_operation_on_iterator<
    Operation, It,
    std::void_t<std::disjunction<is_iterable<It>, std::is_pointer<It>>,
                is_unary_operation<Operation, decltype(*std::declval<It>())>>> = true;

/// Checks whether a string starts with the given prefix
[[nodiscard]] constexpr bool starts_with(std::string_view s, std::string_view prefix) {
  if (s.size() < prefix.size()) {
    return false;
  }
  for (std::size_t i = 0; i < prefix.size(); ++i) {
    if (s[i] != prefix[i]) {
      return false;
    }
  }
  return true;
}

/// Checks whether a string ends with the given prefix
[[nodiscard]] constexpr bool ends_with(std::string_view s, std::string_view suffix) {
  if (s.size() < suffix.size()) {
    return false;
  }
  for (std::size_t i = 1; i <= suffix.size(); ++i) {
    const auto i1 = s.size() - i;
    const auto i2 = suffix.size() - i;
    if (s[i1] != suffix[i2]) {
      return false;
    }
  }
  return true;
}

/// Checks whether a fmt format string starts with the given prefix
template <typename FormatParseContext>
[[nodiscard]] constexpr bool starts_with(const FormatParseContext &ctx, std::string_view prefix) {
  std::string_view format_str{&(*ctx.begin()), static_cast<std::size_t>(ctx.end() - ctx.begin())};
  return starts_with(format_str, prefix);
}

}  // namespace coolerpp
