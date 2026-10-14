// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <stdexcept>
#include <string_view>
#include <utility>

#include "hictk/bin.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/hash.hpp"
#include "hictk/numeric_utils.hpp"

// NOLINTNEXTLINE(*-concat-nested-namespaces)
namespace hictk {

template <typename N>
ThinPixel<N>::operator bool() const noexcept {
  return bin1_id != null_id && bin2_id != null_id;
}

template <typename N>
bool ThinPixel<N>::operator==(const ThinPixel &other) const noexcept {
  return bin1_id == other.bin1_id && bin2_id == other.bin2_id && count == other.count;
}

template <typename N>
bool ThinPixel<N>::operator!=(const ThinPixel &other) const noexcept {
  return !(*this == other);
}

template <typename N>
bool ThinPixel<N>::operator<(const ThinPixel &other) const noexcept {
  if (bin1_id != other.bin1_id) {
    return bin1_id < other.bin1_id;
  }
  if (bin2_id != other.bin2_id) {
    return bin2_id < other.bin2_id;
  }
  return count < other.count;
}

template <typename N>
bool ThinPixel<N>::operator<=(const ThinPixel &other) const noexcept {
  if (bin1_id != other.bin1_id) {
    return bin1_id <= other.bin1_id;
  }
  if (bin2_id != other.bin2_id) {
    return bin2_id <= other.bin2_id;
  }
  return count <= other.count;
}

template <typename N>
bool ThinPixel<N>::operator>(const ThinPixel &other) const noexcept {
  if (bin1_id != other.bin1_id) {
    return bin1_id > other.bin1_id;
  }
  if (bin2_id != other.bin2_id) {
    return bin2_id > other.bin2_id;
  }
  return count > other.count;
}

template <typename N>
bool ThinPixel<N>::operator>=(const ThinPixel &other) const noexcept {
  if (bin1_id != other.bin1_id) {
    return bin1_id >= other.bin1_id;
  }
  if (bin2_id != other.bin2_id) {
    return bin2_id >= other.bin2_id;
  }
  return count >= other.count;
}

template <typename N>
auto ThinPixel<N>::from_coo(const BinTable &bins, std::string_view line, std::int64_t offset)
    -> ThinPixel {
  try {
    auto tp = ThinPixel::from_coo(line, offset);
    if (tp.bin1_id > bins.size()) {
      throw std::out_of_range("invalid bin1_id: out of range");
    }
    if (tp.bin2_id > bins.size()) {
      throw std::out_of_range("invalid bin2_id: out of range");
    }
    return tp;
  } catch (const std::exception &e) {
    internal::raise_pixel_parsing_error(line, "coo", e.what());
  }
}

template <typename N>
auto ThinPixel<N>::from_coo(std::string_view line, std::int64_t offset) -> ThinPixel {
  try {
    const auto tp = internal::ThinPixel::from_coo(line, offset);
    const auto count = internal::parse_numeric_or_throw<N>(tp.count);

    return {tp.bin1_id, tp.bin2_id, count};
  } catch (const std::exception &e) {
    internal::raise_pixel_parsing_error(line, "coo", e.what());
  }
}

template <typename N>
Pixel<N>::Pixel(const Bin &bin, N count_) noexcept : Pixel(bin, bin, count_) {}

template <typename N>
Pixel<N>::Pixel(Bin bin1_, Bin bin2_, N count_) noexcept
    : Pixel({std::move(bin1_), std::move(bin2_)}, count_) {}

template <typename N>
Pixel<N>::Pixel(PixelCoordinates coords_, N count_) noexcept
    : coords(std::move(coords_)), count(count_) {}

template <typename N>
Pixel<N>::Pixel(const Chromosome &chrom, std::uint32_t start, std::uint32_t end, N count_) noexcept
    : Pixel(chrom, start, end, chrom, start, end, count_) {}

template <typename N>
Pixel<N>::Pixel(const Chromosome &chrom1, std::uint32_t start1, std::uint32_t end1,
                const Chromosome &chrom2, std::uint32_t start2, std::uint32_t end2,
                N count_) noexcept
    : Pixel(Bin{chrom1, start1, end1}, Bin{chrom2, start2, end2}, count_) {}

template <typename N>
Pixel<N>::Pixel(const BinTable &bins, std::uint64_t bin_id, N count_)
    : Pixel(bins.at(bin_id), count_) {}

template <typename N>
Pixel<N>::Pixel(const BinTable &bins, const ThinPixel<N> &p)
    : Pixel(bins, p.bin1_id, p.bin2_id, p.count) {}

template <typename N>
Pixel<N>::Pixel(const BinTable &bins, std::uint64_t bin1_id, std::uint64_t bin2_id, N count_)
    : Pixel(bins.at(bin1_id), bins.at(bin2_id), count_) {}

template <typename N>
Pixel<N>::operator bool() const noexcept {
  return !!coords;
}
template <typename N>
bool Pixel<N>::operator==(const Pixel &other) const noexcept {
  return coords == other.coords && count == other.count;
}
template <typename N>
bool Pixel<N>::operator!=(const Pixel &other) const noexcept {
  return !(*this == other);
}
template <typename N>
bool Pixel<N>::operator<(const Pixel &other) const noexcept {
  if (coords == other.coords) {
    return count < other.count;
  }
  return coords < other.coords;
}
template <typename N>
bool Pixel<N>::operator<=(const Pixel &other) const noexcept {
  if (coords == other.coords) {
    return count <= other.count;
  }
  return coords <= other.coords;
}
template <typename N>
bool Pixel<N>::operator>(const Pixel &other) const noexcept {
  if (coords == other.coords) {
    return count > other.count;
  }
  return coords > other.coords;
}
template <typename N>
bool Pixel<N>::operator>=(const Pixel &other) const noexcept {
  if (coords == other.coords) {
    return count >= other.count;
  }
  return coords >= other.coords;
}

template <typename N>
ThinPixel<N> Pixel<N>::to_thin() const noexcept {
  return ThinPixel<N>{coords.bin1.id(), coords.bin2.id(), count};
}

template <typename N>
auto Pixel<N>::from_coo(const BinTable &bins, std::string_view line, std::int64_t offset) -> Pixel {
  try {
    const auto tp = ThinPixel<N>::from_coo(bins, line, offset);
    return {bins.at(tp.bin1_id), bins.at(tp.bin2_id), tp.count};
  } catch (const std::exception &e) {
    internal::raise_pixel_parsing_error(line, "coo", e.what());
  }
}

template <typename N>
auto Pixel<N>::from_bg2(const BinTable &bins, std::string_view line, std::int64_t offset) -> Pixel {
  try {
    const auto p = internal::Pixel::from_bg2(bins, line, offset);
    const auto count = internal::parse_numeric_or_throw<N>(p.count);
    return {std::move(p.bin1), std::move(p.bin2), count};
  } catch (const std::exception &e) {
    internal::raise_pixel_parsing_error(line, "bedgraph2", e.what());
  }
}

template <typename N>
auto Pixel<N>::from_validpair(const BinTable &bins, std::string_view line, std::int64_t offset)
    -> Pixel {
  try {
    const auto p = internal::Pixel::from_validpair(bins, line, offset);
    assert(p.count.empty());
    return {std::move(p.bin1), std::move(p.bin2), 1};
  } catch (const std::exception &e) {
    internal::raise_pixel_parsing_error(line, "validpair", e.what());
  }
}

template <typename N>
auto Pixel<N>::from_4dn_pairs(const BinTable &bins, std::string_view line, std::int64_t offset)
    -> Pixel {
  try {
    const auto p = internal::Pixel::from_4dn_pairs(bins, line, offset);
    assert(p.count.empty());
    return {std::move(p.bin1), std::move(p.bin2), 1};
  } catch (const std::exception &e) {
    internal::raise_pixel_parsing_error(line, "4DN-DCIC pair", e.what());
  }
}

}  // namespace hictk

template <typename N>
std::size_t std::hash<hictk::ThinPixel<N>>::operator()(
    const hictk::ThinPixel<N> &p) const noexcept {
  return hictk::internal::hash_combine(0, p.bin1_id, p.bin2_id, p.count);
}

template <typename N>
std::size_t std::hash<hictk::Pixel<N>>::operator()(const hictk::Pixel<N> &p) const noexcept {
  return hictk::internal::hash_combine(0, p.coords.bin1.id(), p.coords.bin2.id(), p.count);
}
