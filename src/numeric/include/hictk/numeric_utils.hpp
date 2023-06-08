// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>
#include <system_error>

namespace hictk::internal {

template <typename N>
void throw_except_from_errc(std::string_view tok, std::size_t idx, const N &field, const char *c,
                            std::errc e);

auto from_chars(const char *first, const char *last, long double &value) noexcept;

template <typename N>
auto from_chars(const char *first, const char *last, N &value) noexcept;

template <typename N>
void parse_numeric_or_throw(std::string_view tok, N &field);

template <typename N>
N parse_numeric_or_throw(std::string_view tok);
}  // namespace hictk::internal

#include "../../numeric_utils_impl.hpp"
