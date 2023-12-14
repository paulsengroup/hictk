// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <limits>
#include <string_view>
#include <type_traits>
#include <utility>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"

namespace hictk {

class Chromosome;

template <typename N>
struct ThinPixel {
  static constexpr auto null_id = std::numeric_limits<std::uint64_t>::max();
  std::uint64_t bin1_id{null_id};  // NOLINT
  std::uint64_t bin2_id{null_id};  // NOLINT
  N count{};                       // NOLINT

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool operator==(const ThinPixel &other) const noexcept;
  [[nodiscard]] bool operator!=(const ThinPixel &other) const noexcept;
  [[nodiscard]] bool operator<(const ThinPixel &other) const noexcept;
  [[nodiscard]] bool operator<=(const ThinPixel &other) const noexcept;
  [[nodiscard]] bool operator>(const ThinPixel &other) const noexcept;
  [[nodiscard]] bool operator>=(const ThinPixel &other) const noexcept;

  static auto from_coo(std::string_view line) -> ThinPixel;
  static auto from_coo(const BinTable &bins, std::string_view line) -> ThinPixel;
};

struct PixelCoordinates {
  Bin bin1;  // NOLINT
  Bin bin2;  // NOLINT

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

  [[nodiscard]] bool is_intra() const noexcept;
};

template <typename N>
struct Pixel {
  static_assert(std::is_arithmetic_v<N>);

  PixelCoordinates coords{};  // NOLINT
  N count{};                  // NOLINT

  Pixel() = default;
  explicit Pixel(Bin bin, N count_ = 0) noexcept;
  Pixel(Bin bin1_, Bin bin2_, N count_ = 0) noexcept;
  explicit Pixel(PixelCoordinates coords_, N count_ = 0) noexcept;
  Pixel(const Chromosome &chrom, std::uint32_t start, std::uint32_t end, N count_ = 0) noexcept;
  Pixel(const Chromosome &chrom1, std::uint32_t start1, std::uint32_t end1,
        const Chromosome &chrom2, std::uint32_t start2, std::uint32_t end2, N count_ = 0) noexcept;

  Pixel(const BinTable &bins, std::uint64_t bin1_id, std::uint64_t bin2_id, N count_ = 0);
  Pixel(const BinTable &bins, std::uint64_t bin_id, N count_ = 0);
  Pixel(const BinTable &bins, const ThinPixel<N> &p);

  [[nodiscard]] explicit operator bool() const noexcept;
  [[nodiscard]] bool operator==(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator!=(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator<(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator<=(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator>(const Pixel<N> &other) const noexcept;
  [[nodiscard]] bool operator>=(const Pixel<N> &other) const noexcept;

  [[nodiscard]] ThinPixel<N> to_thin() const noexcept;
  static auto from_coo(const BinTable &bins, std::string_view line) -> Pixel;
  static auto from_bg2(const BinTable &bins, std::string_view line) -> Pixel;
  static auto from_validpair(const BinTable &bins, std::string_view line) -> Pixel;
  static auto from_4dn_pairs(const BinTable &bins, std::string_view line) -> Pixel;
};

}  // namespace hictk

#include "./impl/pixel_impl.hpp"  // NOLINT
