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
  enum class Presentation : std::uint_fast8_t { bed, ucsc };
  Presentation presentation{Presentation::ucsc};

  constexpr format_parse_context::iterator parse(format_parse_context &ctx) {
    const auto *it = ctx.begin();
    const auto *end = ctx.end();

    if (hictk::internal::starts_with(ctx, "bed")) {
      presentation = Presentation::bed;
      it += std::string_view{"bed"}.size();  // NOLINT
    } else if (hictk::internal::starts_with(ctx, "ucsc")) {
      presentation = Presentation::ucsc;
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

    if (presentation == Presentation::bed) {
      return format_to(ctx.out(), FMT_STRING("{}\t{}\t{}"), name, gi.start(), gi.end());
    }
    assert(presentation == Presentation::ucsc);
    return format_to(ctx.out(), FMT_STRING("{}:{}-{}"), name, gi.start(), gi.end());
  }
};
}  // namespace fmt
