// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include "hictk/fmt/common.hpp"
#include "hictk/fmt/genomic_interval.hpp"

namespace fmt {
template <>
struct formatter<hictk::Bin> {
  enum Presentation { bed, raw, ucsc };
  Presentation presentation{Presentation::raw};

  constexpr format_parse_context::iterator parse(format_parse_context &ctx) {
    const auto *it = ctx.begin();
    const auto *end = ctx.end();

    if (hictk::internal::starts_with(ctx, "bed")) {
      this->presentation = Presentation::bed;
      it += std::string_view{"bed"}.size();  // NOLINT
    } else if (hictk::internal::starts_with(ctx, "raw")) {
      this->presentation = Presentation::raw;
      it += std::string_view{"raw"}.size();  // NOLINT
    } else if (hictk::internal::starts_with(ctx, "ucsc")) {
      this->presentation = Presentation::ucsc;
      it += std::string_view{"ucsc"}.size();  // NOLINT
    }

    if (it != end && *it != '}') {
      throw fmt::format_error("invalid format");
    }

    return it;
  }
  format_context::iterator format(const hictk::Bin &b, format_context &ctx) const {
    if (this->presentation == Presentation::bed) {
      return fmt::format_to(ctx.out(), FMT_STRING("{:bed}"), b.interval());
    }
    if (this->presentation == Presentation::raw) {
      return fmt::format_to(ctx.out(), FMT_STRING("{}"), b.id());
    }
    assert(this->presentation == Presentation::ucsc);
    return fmt::format_to(ctx.out(), FMT_STRING("{:ucsc}"), b.interval());
  }
};
}  // namespace fmt
