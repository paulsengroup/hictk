// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <string_view>

#include "./common.hpp"
#include "hictk/genomic_interval.hpp"

namespace fmt {
template <>
struct formatter<hictk::GenomicInterval> {
  enum Presentation : std::uint_fast8_t { bed, ucsc };
  Presentation presentation{ucsc};

  constexpr format_parse_context::iterator parse(format_parse_context &ctx) {
    const auto *it = ctx.begin();
    const auto *end = ctx.end();

    if (hictk::internal::starts_with(ctx, "bed")) {
      presentation = bed;
      it += std::string_view{"bed"}.size();  // NOLINT
    } else if (hictk::internal::starts_with(ctx, "ucsc")) {
      presentation = ucsc;
      it += std::string_view{"ucsc"}.size();  // NOLINT
    }

    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    return it;
  }

  template <typename FormatContext>
  auto format(const hictk::GenomicInterval &gi, FormatContext &ctx) const {
    const std::string_view name = !gi ? "null" : gi.chrom().name();

    if (presentation == bed) {
      return format_to(ctx.out(), FMT_STRING("{}\t{}\t{}"), name, gi.start(), gi.end());
    }
    assert(presentation == ucsc);
    return format_to(ctx.out(), FMT_STRING("{}:{}-{}"), name, gi.start(), gi.end());
  }
};
}  // namespace fmt
