// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <string_view>

#include "hictk/fmt/bin.hpp"
#include "hictk/fmt/common.hpp"
#include "hictk/pixel.hpp"

template <>
struct fmt::formatter<hictk::PixelCoordinates> {
  enum class Presentation : std::uint_fast8_t { bg2, raw };
  Presentation presentation{Presentation::bg2};

  constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    const auto *it = ctx.begin();
    const auto *end = ctx.end();

    if (hictk::internal::starts_with(ctx, "bg2")) {
      presentation = Presentation::bg2;
      it += std::string_view{"bg2"}.size();  // NOLINT
    } else if (hictk::internal::starts_with(ctx, "raw")) {
      presentation = Presentation::raw;
      it += std::string_view{"raw"}.size();  // NOLINT
    }

    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    return it;
  }

  template <typename FormatContext>
  auto format(const hictk::PixelCoordinates &c, FormatContext &ctx) const -> decltype(ctx.out()) {
    if (presentation == Presentation::bg2) {
      return fmt::format_to(ctx.out(), FMT_STRING("{:bed}\t{:bed}"), c.bin1, c.bin2);
    }

    assert(presentation == Presentation::raw);
    return fmt::format_to(ctx.out(), FMT_STRING("{:raw}\t{:raw}"), c.bin1, c.bin2);
  }
};

template <typename N>
struct fmt::formatter<hictk::Pixel<N>> {
  enum class Presentation : std::uint_fast8_t { bg2, raw };
  Presentation presentation{Presentation::bg2};

  constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
    const auto *it = ctx.begin();
    const auto *end = ctx.end();

    if (hictk::internal::starts_with(ctx, "bg2")) {
      presentation = Presentation::bg2;
      it += std::string_view{"bg2"}.size();  // NOLINT
    } else if (hictk::internal::starts_with(ctx, "raw")) {
      presentation = Presentation::raw;
      it += std::string_view{"raw"}.size();  // NOLINT
    }

    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    return it;
  }

  template <typename FormatContext>
  auto format(const hictk::Pixel<N> &p, FormatContext &ctx) const -> decltype(ctx.out()) {
    if (presentation == Presentation::raw) {
      return fmt::format_to(ctx.out(), FMT_STRING("{:raw}\t{}"), p.coords, p.count);
    }

    assert(presentation == Presentation::bg2);
    return fmt::format_to(ctx.out(), FMT_STRING("{:bg2}\t{}"), p.coords, p.count);
  }
};

template <typename N>
struct fmt::formatter<hictk::ThinPixel<N>> {
  static constexpr auto parse(const format_parse_context &ctx) -> decltype(ctx.begin()) {
    const auto *it = ctx.begin();
    const auto *end = ctx.end();
    if (it != end && *it != '}') {
      throw format_error("invalid format");
    }

    return it;
  }

  template <typename FormatContext>
  auto format(const hictk::ThinPixel<N> &p, FormatContext &ctx) const -> decltype(ctx.out()) {
    return fmt::format_to(ctx.out(), FMT_STRING("{:d}\t{:d}\t{}"), p.bin1_id, p.bin2_id, p.count);
  }
};
