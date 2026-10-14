// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/pixel.hpp"

#include <fmt/format.h>

#include <array>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string_view>
#include <utility>

#include "hictk/bin.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/numeric_utils.hpp"

namespace hictk {

namespace internal {

template <std::size_t N>
[[nodiscard]] static std::array<std::string_view, N> tokenize_n(std::string_view line) {
  std::array<std::string_view, N> toks{};
  for (std::size_t i = 0; i < toks.size(); ++i) {
    const auto pos = line.find('\t');
    toks[i] = line.substr(0, pos);
    if (toks[i].empty()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("expected exactly {} fields, found {}"), N, i));
    }

    if (pos == std::string_view::npos) {
      line = "";
    } else {
      line.remove_prefix(pos + 1);
    }
  }

  if (!line.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("expected exactly {} fields, found {} or more"), N, N + 1));
  }

  return toks;
}

template <std::size_t N>
[[nodiscard]] static std::array<std::string_view, N> tokenize_n_or_more(std::string_view line) {
  std::array<std::string_view, N> toks{};
  for (std::size_t i = 0; i < toks.size(); ++i) {
    const auto pos = line.find('\t');
    toks[i] = line.substr(0, pos);
    if (toks[i].empty()) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("expected {} or more fields, found {}"), N, i));
    }

    if (pos == std::string_view::npos) {
      line = "";
    } else {
      line.remove_prefix(pos + 1);
    }
  }
  return toks;
}

auto ThinPixel::from_coo(std::string_view line, std::int64_t offset) -> ThinPixel {
  if (line.empty()) {
    throw std::runtime_error("found an empty line");
  }
  if (line.back() == '\r') {
    line = line.substr(0, line.size() - 1);
  }
  const auto toks = internal::tokenize_n<3>(line);

  const auto bin1_id =
      static_cast<std::uint64_t>(internal::parse_numeric_or_throw<std::int64_t>(toks[0]) + offset);
  const auto bin2_id =
      static_cast<std::uint64_t>(internal::parse_numeric_or_throw<std::int64_t>(toks[1]) + offset);
  const auto &count = toks[2];  // NOLINT(*-avoid-magic-numbers)

  return {bin1_id, bin2_id, count};
}

auto Pixel::from_bg2(const BinTable &bins, std::string_view line, std::int64_t offset) -> Pixel {
  if (line.empty()) {
    throw std::runtime_error("found an empty line");
  }
  if (line.back() == '\r') {
    line = line.substr(0, line.size() - 1);
  }
  const auto toks = internal::tokenize_n_or_more<7>(line);

  const auto &chrom1 = toks[0];
  const auto start1 =
      static_cast<std::uint32_t>(internal::parse_numeric_or_throw<std::int64_t>(toks[1]) + offset);

  const auto &chrom2 = toks[3];
  const auto start2 =
      static_cast<std::uint32_t>(internal::parse_numeric_or_throw<std::int64_t>(toks[4]) + offset);

  const auto &count = toks[6];  // NOLINT(*-avoid-magic-numbers)
  return {bins.at(chrom1, start1), bins.at(chrom2, start2), count};
}

auto Pixel::from_validpair(const BinTable &bins, std::string_view line, std::int64_t offset)
    -> Pixel {
  if (line.empty()) {
    throw std::runtime_error("found an empty line");
  }
  if (line.back() == '\r') {
    line = line.substr(0, line.size() - 1);
  }
  const auto toks = internal::tokenize_n_or_more<6>(line);

  const auto &chrom1 = toks[1];
  const auto start1 =
      static_cast<std::uint32_t>(internal::parse_numeric_or_throw<std::int64_t>(toks[2]) + offset);

  const auto &chrom2 = toks[4];
  const auto start2 =
      static_cast<std::uint32_t>(internal::parse_numeric_or_throw<std::int64_t>(toks[5]) + offset);

  return {bins.at(chrom1, start1), bins.at(chrom2, start2), ""};
}

auto Pixel::from_4dn_pairs(const BinTable &bins, std::string_view line, std::int64_t offset)
    -> Pixel {
  if (line.empty()) {
    throw std::runtime_error("found an empty line");
  }
  if (line.back() == '\r') {
    line = line.substr(0, line.size() - 1);
  }
  const auto toks = internal::tokenize_n_or_more<6>(line);

  const auto &chrom1 = toks[1];
  const auto start1 =
      static_cast<std::uint32_t>(internal::parse_numeric_or_throw<std::int64_t>(toks[2]) + offset);

  const auto &chrom2 = toks[3];
  const auto start2 =
      static_cast<std::uint32_t>(internal::parse_numeric_or_throw<std::int64_t>(toks[4]) + offset);

  return {bins.at(chrom1, start1), bins.at(chrom2, start2), ""};
}

void raise_pixel_parsing_error(std::string_view line, std::string_view format,
                               std::string_view error_msg) {
  throw std::runtime_error(
      fmt::format(FMT_STRING("line \"{}\" is not in {} format: {}"), line, format, error_msg));
}

}  // namespace internal

PixelCoordinates::PixelCoordinates(Bin bin1_, Bin bin2_) noexcept
    : bin1(std::move(bin1_)), bin2(std::move(bin2_)) {}

PixelCoordinates::PixelCoordinates(std::pair<Bin, Bin> bins) noexcept
    : PixelCoordinates(std::move(bins.first), std::move(bins.second)) {}

PixelCoordinates::PixelCoordinates(Bin bin) noexcept : bin1(bin), bin2(std::move(bin)) {}

PixelCoordinates::operator bool() const noexcept { return !empty(); }

bool PixelCoordinates::empty() const noexcept { return !bin1 && !bin2; }

bool PixelCoordinates::operator==(const PixelCoordinates &other) const noexcept {
  return bin1 == other.bin1 && bin2 == other.bin2;
}

bool PixelCoordinates::operator!=(const PixelCoordinates &other) const noexcept {
  return !(*this == other);
}

bool PixelCoordinates::operator<(const PixelCoordinates &other) const noexcept {
  if (bin1 == other.bin1) {
    return bin2 < other.bin2;
  }
  return bin1 < other.bin1;
}

bool PixelCoordinates::operator<=(const PixelCoordinates &other) const noexcept {
  if (bin1 == other.bin1) {
    return bin2 <= other.bin2;
  }
  return bin1 <= other.bin1;
}

bool PixelCoordinates::operator>(const PixelCoordinates &other) const noexcept {
  if (bin1 == other.bin1) {
    return bin2 > other.bin2;
  }
  return bin1 > other.bin1;
}

bool PixelCoordinates::operator>=(const PixelCoordinates &other) const noexcept {
  if (bin1 == other.bin1) {
    return bin2 >= other.bin2;
  }
  return bin1 >= other.bin1;
}

bool PixelCoordinates::is_intra() const noexcept { return bin1.chrom() == bin2.chrom(); }

std::uint64_t area(const PixelCoordinates &coords, std::uint32_t resolution,
                   bool upper_triangular) noexcept {
  return area(coords, coords, resolution, upper_triangular);
}

std::uint64_t area(const PixelCoordinates &coords1, const PixelCoordinates &coords2,
                   std::uint32_t resolution, bool upper_triangular) noexcept {
  return area(GenomicInterval{coords1.bin1.chrom(), coords1.bin1.start(), coords1.bin2.end()},
              GenomicInterval{coords2.bin1.chrom(), coords2.bin1.start(), coords2.bin2.end()},
              resolution, upper_triangular);
}

}  // namespace hictk
