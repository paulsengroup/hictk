// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <string_view>

#include "hictk/bin.hpp"
#include "hictk/fmt/common.hpp"
#include "hictk/fmt/genomic_interval.hpp"

namespace fmt {
template <>
struct formatter<hictk::Bin> {
  enum Presentation : std::uint_fast8_t { bed, raw, ucsc };
  Presentation presentation{raw};

  constexpr format_parse_context::iterator parse(format_parse_context &ctx) {
    const auto *it = ctx.begin();
    const auto *end = ctx.end();

    // NOLINTBEGIN(*-bounds-pointer-arithmetic)
    if (hictk::internal::starts_with(ctx, "bed")) {
      presentation = bed;
      it += std::string_view{"bed"}.size();  // NOLINT
    } else if (hictk::internal::starts_with(ctx, "raw")) {
      presentation = raw;
      it += std::string_view{"raw"}.size();  // NOLINT
    } else if (hictk::internal::starts_with(ctx, "ucsc")) {
      presentation = ucsc;
      it += std::string_view{"ucsc"}.size();  // NOLINT
    }
    // NOLINTEND(*-bounds-pointer-arithmetic)

    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    return it;
  }
  format_context::iterator format(const hictk::Bin &b, format_context &ctx) const {
    if (presentation == bed) {
      return format_to(ctx.out(), FMT_STRING("{:bed}"), b.interval());
    }
    if (presentation == raw) {
      return format_to(ctx.out(), FMT_STRING("{}"), b.id());
    }
    assert(presentation == ucsc);
    return format_to(ctx.out(), FMT_STRING("{:ucsc}"), b.interval());
  }
};
}  // namespace fmt
