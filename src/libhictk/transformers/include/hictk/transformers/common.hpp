// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <type_traits>

namespace hictk::transformers {
enum class QuerySpan : std::uint_fast8_t { lower_triangle, upper_triangle, full };

namespace internal {

template <typename T, typename = std::void_t<>>
inline constexpr bool has_coord1_member_fx = false;

template <typename T>
inline constexpr bool has_coord1_member_fx<T, std::void_t<decltype(std::declval<T>().coord1())>> =
    true;

}  // namespace internal
}  // namespace hictk::transformers
