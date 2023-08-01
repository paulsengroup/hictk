// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include "hictk/fmt/common.hpp"
#include "hictk/genomic_interval.hpp"

namespace fmt {
template <>
struct formatter<hictk::GenomicInterval> {
  enum Presentation { bed, ucsc };
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
      throw fmt::format_error("invalid format");
    }

    return it;
  }
  inline format_context::iterator format(const hictk::GenomicInterval &gi,
                                         format_context &ctx) const {
    const std::string_view name = !gi ? "null" : gi.chrom().name();

    if (presentation == Presentation::bed) {
      return fmt::format_to(ctx.out(), FMT_STRING("{}\t{}\t{}"), name, gi.start(), gi.end());
    }
    assert(presentation == Presentation::ucsc);
    return fmt::format_to(ctx.out(), FMT_STRING("{}:{}-{}"), name, gi.start(), gi.end());
  }
};
}  // namespace fmt
