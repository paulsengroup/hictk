// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <string_view>

#include "./common.hpp"

namespace fmt {
template <>
struct formatter<hictk::Chromosome> {
  enum class Presentation : std::uint_fast8_t { tsv, ucsc };
  Presentation presentation{Presentation::ucsc};

  constexpr format_parse_context::iterator parse(format_parse_context& ctx) {
    const auto* it = ctx.begin();
    const auto* end = ctx.end();

    // NOLINTBEGIN(*-bounds-pointer-arithmetic)
    if (hictk::internal::starts_with(ctx, "ucsc")) {
      presentation = Presentation::ucsc;
      it += std::string_view{"ucsc"}.size();
    } else if (hictk::internal::starts_with(ctx, "tsv")) {
      presentation = Presentation::tsv;
      it += std::string_view{"tsv"}.size();
    }
    // NOLINTEND(*-bounds-pointer-arithmetic)

    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    return it;
  }

  template <typename FormatContext>
  auto format(const hictk::Chromosome& c, FormatContext& ctx) const {
    return presentation == Presentation::tsv
               ? format_to(ctx.out(), FMT_STRING("{}\t{}"), c.name(), c.size())
               : format_to(ctx.out(), FMT_STRING("{}:{}"), c.name(), c.size());
  }
};  // namespace fmt

}  // namespace fmt
