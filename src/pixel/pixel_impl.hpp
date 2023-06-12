// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <cstdint>
#include <string_view>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"

namespace hictk {
/*
inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTable> &bins,
                                          const Chromosome &chrom1, const Chromosome &chrom2,
                                          std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->chromosomes().get_id(chrom1), bins->chromosomes().get_id(chrom2),
                       bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTable> &bins,
                                          std::string_view chrom1_name,
                                          std::string_view chrom2_name, std::uint32_t bin1_start_,
                                          std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->chromosomes().get_id(chrom1_name),
                       bins->chromosomes().get_id(chrom2_name), bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTable> &bins,
                                          std::uint32_t chrom1_id_, std::uint32_t chrom2_id_,
                                          std::uint32_t bin1_start_, std::uint32_t bin2_start_)
    : PixelCoordinates(bins, bins->map_to_bin_id(chrom1_id_, bin1_start_),
                       bins->map_to_bin_id(chrom2_id_, bin2_start_)) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTable> &bins,
                                          const Chromosome &chrom, std::uint32_t bin1_start_,
                                          std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom, chrom, bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTable> &bins,
                                          std::uint32_t chrom_id, std::uint32_t bin1_start_,
                                          std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom_id, chrom_id, bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(const std::shared_ptr<const BinTable> &bins,
                                          std::string_view chrom_name, std::uint32_t bin1_start_,
                                          std::uint32_t bin2_start_)
    : PixelCoordinates(bins, chrom_name, chrom_name, bin1_start_, bin2_start_) {}

inline PixelCoordinates::PixelCoordinates(std::shared_ptr<const BinTable> bins,
                                          std::uint64_t bin1_id_, std::uint64_t bin2_id_)
    : _bins(std::move(bins)), _bin1_id(bin1_id_), _bin2_id(bin2_id_) {
  assert(_bin1_id <= _bins->size());
  assert(_bin2_id <= _bins->size());
}

inline PixelCoordinates::operator bool() const noexcept { return !!this->_bins; }

inline const Chromosome &PixelCoordinates::chrom1() const { return this->bin1().chrom; }

inline const Chromosome &PixelCoordinates::chrom2() const { return this->bin2().chrom; }

inline std::uint32_t PixelCoordinates::bin1.chrom().id() const {
  return this->_bins->chromosomes().get_id(this->chrom1());
}

inline std::uint32_t PixelCoordinates::bin2.chrom().id() const {
  return this->_bins->chromosomes().get_id(this->chrom2());
}

inline GenomicInterval PixelCoordinates::bin1() const {
  assert(this->_bins);
  assert(!!*this);

  return this->_bins->bin_id_to_coords(_bin1_id);
}

inline GenomicInterval PixelCoordinates::bin2() const {
  assert(this->_bins);
  assert(!!*this);

  return this->_bins->at(_bin2_id);
}

inline std::uint64_t PixelCoordinates::bin1_id() const noexcept { return this->_bin1_id; }
inline std::uint64_t PixelCoordinates::bin2_id() const noexcept { return this->_bin2_id; }

inline std::uint32_t PixelCoordinates::bin_size() const noexcept {
  assert(this->_bins);
  return this->_bins->bin_size();
}

constexpr bool PixelCoordinates::operator==(const PixelCoordinates &other) const noexcept {
  return this->_bin1_id == other._bin1_id && this->_bin2_id == other._bin2_id;
}

constexpr bool PixelCoordinates::operator!=(const PixelCoordinates &other) const noexcept {
  return !(*this == other);
}

constexpr bool PixelCoordinates::operator<(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id < other._bin2_id;
  }
  return this->_bin1_id < other._bin1_id;
}

constexpr bool PixelCoordinates::operator<=(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id <= other._bin2_id;
  }
  return this->_bin1_id <= other._bin1_id;
}

constexpr bool PixelCoordinates::operator>(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id > other._bin2_id;
  }
  return this->_bin1_id > other._bin1_id;
}

constexpr bool PixelCoordinates::operator>=(const PixelCoordinates &other) const noexcept {
  if (this->_bin1_id == other._bin1_id) {
    return this->_bin2_id >= other._bin2_id;
  }
  return this->_bin1_id >= other._bin1_id;
}
*/

inline PixelCoordinates::PixelCoordinates(Bin bin1_, Bin bin2_) noexcept
    : bin1(std::move(bin1_)), bin2(std::move(bin2_)) {}

inline PixelCoordinates::PixelCoordinates(std::pair<Bin, Bin> bins) noexcept
    : PixelCoordinates(std::move(bins.first), std::move(bins.second)) {}

inline PixelCoordinates::PixelCoordinates(Bin bin) noexcept : bin1(bin), bin2(std::move(bin)) {}

inline PixelCoordinates::operator bool() const noexcept { return !!this->bin1 && !!this->bin2; }

inline bool PixelCoordinates::operator==(const PixelCoordinates &other) const noexcept {
  return this->bin1 == other.bin1 && this->bin2 == other.bin2;
}

inline bool PixelCoordinates::operator!=(const PixelCoordinates &other) const noexcept {
  return !(*this == other);
}

inline bool PixelCoordinates::operator<(const PixelCoordinates &other) const noexcept {
  if (this->bin1 == other.bin1) {
    return this->bin2 < other.bin2;
  }
  return this->bin1 < other.bin1;
}

inline bool PixelCoordinates::operator<=(const PixelCoordinates &other) const noexcept {
  if (this->bin1 == other.bin1) {
    return this->bin2 <= other.bin2;
  }
  return this->bin1 <= other.bin1;
}

inline bool PixelCoordinates::operator>(const PixelCoordinates &other) const noexcept {
  if (this->bin1 == other.bin1) {
    return this->bin2 > other.bin2;
  }
  return this->bin1 > other.bin1;
}

inline bool PixelCoordinates::operator>=(const PixelCoordinates &other) const noexcept {
  if (this->bin1 == other.bin1) {
    return this->bin2 >= other.bin2;
  }
  return this->bin1 >= other.bin1;
}

inline bool PixelCoordinates::is_intra() const noexcept {
  return this->bin1.chrom() == this->bin2.chrom();
}

template <typename N>
inline Pixel<N>::Pixel(Bin bin, N count_) noexcept : Pixel(bin, std::move(bin), count_) {}

template <typename N>
inline Pixel<N>::Pixel(Bin bin1_, Bin bin2_, N count_) noexcept
    : Pixel({std::move(bin1_), std::move(bin2_)}, count_) {}

template <typename N>
inline Pixel<N>::Pixel(PixelCoordinates coords_, N count_) noexcept
    : coords(std::move(coords_)), count(count_) {}

template <typename N>
inline Pixel<N>::Pixel(const Chromosome &chrom, std::uint32_t start, std::uint32_t end,
                       N count_) noexcept
    : Pixel(chrom, start, end, chrom, start, end, count_) {}

template <typename N>
inline Pixel<N>::Pixel(const Chromosome &chrom1, std::uint32_t start1, std::uint32_t end1,
                       const Chromosome &chrom2, std::uint32_t start2, std::uint32_t end2,
                       N count_) noexcept
    : Pixel(Bin{chrom1, start1, end1}, Bin{chrom2, start2, end2}, count_) {}

template <typename N>
inline Pixel<N>::Pixel(const BinTable &bins, std::uint64_t bin_id, N count_)
    : Pixel(bins.at(bin_id), count_) {}

template <typename N>
inline Pixel<N>::Pixel(const BinTable &bins, std::uint64_t bin1_id, std::uint64_t bin2_id, N count_)
    : Pixel(bins.at(bin1_id), bins.at(bin2_id), count_) {}

template <typename N>
inline Pixel<N>::operator bool() const noexcept {
  return !!this->coords;
}
template <typename N>
inline bool Pixel<N>::operator==(const Pixel<N> &other) const noexcept {
  return this->coords == other.coords && this->count == other.count;
}
template <typename N>
inline bool Pixel<N>::operator!=(const Pixel<N> &other) const noexcept {
  return !(*this == other);
}
template <typename N>
inline bool Pixel<N>::operator<(const Pixel<N> &other) const noexcept {
  if (this->coords == other.coords) {
    return this->count < other.count;
  }
  return this->coords < other.coords;
}
template <typename N>
inline bool Pixel<N>::operator<=(const Pixel<N> &other) const noexcept {
  if (this->coords == other.coords) {
    return this->count <= other.count;
  }
  return this->coords <= other.coords;
}
template <typename N>
inline bool Pixel<N>::operator>(const Pixel<N> &other) const noexcept {
  if (this->coords == other.coords) {
    return this->count > other.count;
  }
  return this->coords > other.coords;
}
template <typename N>
inline bool Pixel<N>::operator>=(const Pixel<N> &other) const noexcept {
  if (this->coords == other.coords) {
    return this->count >= other.count;
  }
  return this->coords >= other.coords;
}
}  // namespace hictk
