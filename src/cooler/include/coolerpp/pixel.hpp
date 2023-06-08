// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cstdint>
#include <limits>
#include <string_view>
#include <type_traits>

#include "coolerpp/bin_table.hpp"

namespace coolerpp {

class Chromosome;

struct PixelCoordinates {
  Bin bin1;
  Bin bin2;

  PixelCoordinates() = default;
  PixelCoordinates(Bin bin1_, Bin bin2_) noexcept;
  explicit PixelCoordinates(std::pair<Bin, Bin> bins) noexcept;
  explicit PixelCoordinates(Bin bin) noexcept;

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool operator==(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator!=(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator<(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator<=(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator>(const PixelCoordinates &other) const noexcept;
  [[nodiscard]] bool operator>=(const PixelCoordinates &other) const noexcept;
};

template <typename N>
struct Pixel {
  static_assert(std::is_arithmetic_v<N>);

  PixelCoordinates coords{};
  N count{};

  Pixel() = default;
  explicit Pixel(Bin bin, N count_ = 0) noexcept;
  Pixel(Bin bin1_, Bin bin2_, N count_ = 0) noexcept;
  explicit Pixel(PixelCoordinates coords_, N count_ = 0) noexcept;
  Pixel(const Chromosome &chrom, std::uint32_t start, std::uint32_t end, N count_ = 0) noexcept;
  Pixel(const Chromosome &chrom1, std::uint32_t start1, std::uint32_t end1,
        const Chromosome &chrom2, std::uint32_t start2, std::uint32_t end2, N count_ = 0) noexcept;

  Pixel(const BinTable &bins, std::uint64_t bin1_id, std::uint64_t bin2_id, N count_ = 0);
  Pixel(const BinTable &bins, std::uint64_t bin_id, N count_ = 0);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool operator==(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator!=(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator<(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator<=(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator>(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator>=(const Pixel<N> &other) const noexcept;
};

}  // namespace coolerpp

template <>
struct fmt::formatter<coolerpp::PixelCoordinates> {
  enum Presentation { bg2, raw };
  Presentation presentation{Presentation::bg2};

  constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin());
  template <typename FormatContext>
  auto format(const coolerpp::PixelCoordinates &c, FormatContext &ctx) const -> decltype(ctx.out());
};

template <typename N>
struct fmt::formatter<coolerpp::Pixel<N>> {
  enum Presentation { bg2, raw };
  Presentation presentation{Presentation::bg2};

  constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin());
  template <typename FormatContext>
  auto format(const coolerpp::Pixel<N> &p, FormatContext &ctx) const -> decltype(ctx.out());
};

#include "../../pixel_impl.hpp"
