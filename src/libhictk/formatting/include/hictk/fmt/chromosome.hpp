// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include "hictk/fmt/common.hpp"

namespace fmt {
template <>
struct formatter<hictk::Chromosome> {
  enum Presentation { tsv, ucsc };
  Presentation presentation{Presentation::ucsc};

  constexpr format_parse_context::iterator parse(format_parse_context& ctx) {
    auto* it = ctx.begin();
    const auto* end = ctx.end();

    if (hictk::internal::starts_with(ctx, "ucsc")) {
      presentation = Presentation::ucsc;
      it += std::string_view{"ucsc"}.size();
    } else if (hictk::internal::starts_with(ctx, "tsv")) {
      presentation = Presentation::tsv;
      it += std::string_view{"tsv"}.size();
    }

    if (it != end && *it != '}') {
      throw fmt::format_error("invalid format");
    }

    return it;
  }

  inline format_context::iterator format(const hictk::Chromosome& c, format_context& ctx) const {
    return presentation == Presentation::tsv
               ? fmt::format_to(ctx.out(), FMT_STRING("{}\t{}"), c.name(), c.size())
               : fmt::format_to(ctx.out(), FMT_STRING("{}:{}"), c.name(), c.size());
  }
};  // namespace fmt

}  // namespace fmt
