// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

namespace hictk::cooler::utils {

constexpr ValidationStatusCooler::operator bool() const noexcept { return is_cooler; }

constexpr ValidationStatusMultiresCooler::operator bool() const noexcept {
  return is_multires_file;
}

constexpr ValidationStatusScool::operator bool() const noexcept { return is_scool_file; }

}  // namespace hictk::cooler::utils

constexpr auto fmt::formatter<hictk::cooler::utils::ValidationStatusCooler>::parse(
    format_parse_context &ctx) -> format_parse_context::iterator {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw format_error("invalid format");
  }
  return ctx.end();
}

constexpr auto fmt::formatter<hictk::cooler::utils::ValidationStatusMultiresCooler>::parse(
    format_parse_context &ctx) -> format_parse_context::iterator {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw format_error("invalid format");
  }
  return ctx.end();
}

constexpr auto fmt::formatter<hictk::cooler::utils::ValidationStatusScool>::parse(
    format_parse_context &ctx) -> format_parse_context::iterator {
  if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
    throw format_error("invalid format");
  }
  return ctx.end();
}
