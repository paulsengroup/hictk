// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <type_traits>

namespace hictk::transformers::internal {
template <typename T, typename = std::void_t<>>
inline constexpr bool has_coord1_member_fx = false;

template <typename T>
inline constexpr bool has_coord1_member_fx<T, std::void_t<decltype(std::declval<T>().coord1())>> =
    true;
}
