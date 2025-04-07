// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstddef>
#include <string_view>

#include "hictk/string.hpp"

namespace hictk::internal {

/// Checks whether a fmt format string starts with the given prefix
[[nodiscard]] constexpr bool starts_with(const fmt::format_parse_context &ctx,
                                         std::string_view prefix) {
  const std::string_view format_str{&(*ctx.begin()),
                                    static_cast<std::size_t>(ctx.end() - ctx.begin())};
  return starts_with(format_str, prefix);
}
}  // namespace hictk::internal
