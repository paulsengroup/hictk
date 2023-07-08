// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <cassert>
#include <cstdint>
#include <string_view>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"

namespace hictk {

namespace internal {
template <std::size_t N>
[[nodiscard]] std::array<std::string_view, N> tokenize_n(std::string_view line) {
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
[[nodiscard]] std::array<std::string_view, N> tokenize_n_or_more(std::string_view line) {
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
}  // namespace internal

template <typename N>
inline ThinPixel<N>::operator bool() const noexcept {
  return this->bin1_id != null_id && this->bin2_id != null_id;
}

template <typename N>
inline bool ThinPixel<N>::operator==(const ThinPixel &other) const noexcept {
  return this->bin1_id == other.bin1_id && this->bin2_id == other.bin2_id &&
         this->count == other.count;
}

template <typename N>
inline bool ThinPixel<N>::operator!=(const ThinPixel &other) const noexcept {
  return !(*this == other);
}

template <typename N>
inline bool ThinPixel<N>::operator<(const ThinPixel &other) const noexcept {
  if (this->bin1_id != other.bin1_id) {
    return this->bin1_id < other.bin1_id;
  }
  if (this->bin2_id != other.bin2_id) {
    return this->bin2_id < other.bin2_id;
  }
  return this->count < other.count;
}

template <typename N>
inline bool ThinPixel<N>::operator<=(const ThinPixel &other) const noexcept {
  if (this->bin1_id != other.bin1_id) {
    return this->bin1_id <= other.bin1_id;
  }
  if (this->bin2_id != other.bin2_id) {
    return this->bin2_id <= other.bin2_id;
  }
  return this->count <= other.count;
}

template <typename N>
inline bool ThinPixel<N>::operator>(const ThinPixel &other) const noexcept {
  if (this->bin1_id != other.bin1_id) {
    return this->bin1_id > other.bin1_id;
  }
  if (this->bin2_id != other.bin2_id) {
    return this->bin2_id > other.bin2_id;
  }
  return this->count > other.count;
}

template <typename N>
inline bool ThinPixel<N>::operator>=(const ThinPixel &other) const noexcept {
  if (this->bin1_id != other.bin1_id) {
    return this->bin1_id >= other.bin1_id;
  }
  if (this->bin2_id != other.bin2_id) {
    return this->bin2_id >= other.bin2_id;
  }
  return this->count >= other.count;
}

template <typename N>
inline auto ThinPixel<N>::from_coo(const BinTable &bins, std::string_view line) -> ThinPixel<N> {
  try {
    const auto toks = internal::tokenize_n<3>(line);

    const auto bin1_id = internal::parse_numeric_or_throw<std::size_t>(toks[0]);
    const auto bin2_id = internal::parse_numeric_or_throw<std::size_t>(toks[1]);
    const auto count = internal::parse_numeric_or_throw<N>(toks[2]);

    if (bin1_id > bins.size()) {
      throw std::out_of_range("invalid bin1_id: out of range");
    }
    if (bin2_id > bins.size()) {
      throw std::out_of_range("invalid bin2_id: out of range");
    }

    return {bins.at(bin1_id).id(), bins.at(bin2_id).id(), count};
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("line \"{}\" is not in coo format: {}"), line, e.what()));
  }
}
template <typename N>
inline auto ThinPixel<N>::from_coo(std::string_view line) -> ThinPixel<N> {
  try {
    const auto toks = internal::tokenize_n<3>(line);

    const auto bin1_id = internal::parse_numeric_or_throw<std::size_t>(toks[0]);
    const auto bin2_id = internal::parse_numeric_or_throw<std::size_t>(toks[1]);
    const auto count = internal::parse_numeric_or_throw<N>(toks[2]);

    return {bin1_id, bin2_id, count};
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("line \"{}\" is not in coo format: {}"), line, e.what()));
  }
}

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
inline Pixel<N>::Pixel(const hictk::BinTable &bins, const hictk::ThinPixel<N> &p)
    : Pixel(bins, p.bin1_id, p.bin2_id, p.count) {}

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

template <typename N>
inline ThinPixel<N> Pixel<N>::to_thin() const noexcept {
  return ThinPixel<N>{this->coords.bin1.id(), this->coords.bin2.id(), this->count};
}

template <typename N>
inline auto Pixel<N>::from_coo(const hictk::BinTable &bins, std::string_view line) -> Pixel<N> {
  try {
    const auto toks = internal::tokenize_n<3>(line);

    const auto bin1_id = internal::parse_numeric_or_throw<std::size_t>(toks[0]);
    const auto bin2_id = internal::parse_numeric_or_throw<std::size_t>(toks[1]);
    const auto count = internal::parse_numeric_or_throw<N>(toks[2]);

    return {bins.at(bin1_id), bins.at(bin2_id), count};
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("line \"{}\" is not in coo format: {}"), line, e.what()));
  }
}

template <typename N>
inline auto Pixel<N>::from_bg2(const hictk::BinTable &bins, std::string_view line) -> Pixel<N> {
  try {
    const auto toks = internal::tokenize_n_or_more<7>(line);

    const auto &chrom1 = toks[0];
    const auto start1 = internal::parse_numeric_or_throw<std::uint32_t>(toks[1]);

    const auto &chrom2 = toks[3];
    const auto start2 = internal::parse_numeric_or_throw<std::uint32_t>(toks[4]);

    const auto count = internal::parse_numeric_or_throw<N>(toks[6]);
    return {bins.at(chrom1, start1), bins.at(chrom2, start2), count};
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("line \"{}\" is not in bedgraph2 format: {}"), line, e.what()));
  }
}

template <typename N>
inline auto Pixel<N>::from_validpair(const hictk::BinTable &bins, std::string_view line)
    -> Pixel<N> {
  try {
    const auto toks = internal::tokenize_n_or_more<6>(line);

    const auto &chrom1 = toks[1];
    const auto start1 = internal::parse_numeric_or_throw<std::uint32_t>(toks[2]);

    const auto &chrom2 = toks[4];
    const auto start2 = internal::parse_numeric_or_throw<std::uint32_t>(toks[5]);

    return {bins.at(chrom1, start1), bins.at(chrom2, start2), 1};
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("line \"{}\" is not in validpair format: {}"), line, e.what()));
  }
}

template <typename N>
inline auto Pixel<N>::from_4dn_pairs(const hictk::BinTable &bins, std::string_view line)
    -> Pixel<N> {
  try {
    const auto toks = internal::tokenize_n_or_more<6>(line);

    const auto &chrom1 = toks[1];
    const auto start1 = internal::parse_numeric_or_throw<std::uint32_t>(toks[2]);

    const auto &chrom2 = toks[3];
    const auto start2 = internal::parse_numeric_or_throw<std::uint32_t>(toks[4]);

    return {bins.at(chrom1, start1), bins.at(chrom2, start2), 1};
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("line \"{}\" is not in 4DN-DCIC pair format: {}"), line, e.what()));
  }
}

namespace internal {
template <typename PixelIt>
inline bool PixelMerger<PixelIt>::Node::operator<(const Node &other) const noexcept {
  assert(!!this->pixel);
  assert(!!other.pixel);
  if (this->pixel.bin1_id != other.pixel.bin1_id) {
    return this->pixel.bin1_id < other.pixel.bin1_id;
  }
  return this->pixel.bin2_id < other.pixel.bin2_id;
}

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::Node::operator>(const Node &other) const noexcept {
  assert(!!this->pixel);
  assert(!!other.pixel);
  if (this->pixel.bin1_id != other.pixel.bin1_id) {
    return this->pixel.bin1_id > other.pixel.bin1_id;
  }
  return this->pixel.bin2_id > other.pixel.bin2_id;
}

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::Node::operator==(const Node &other) const noexcept {
  return this->pixel.bin1_id == other.pixel.bin1_id && this->pixel.bin2_id == other.pixel.bin2_id;
}

template <typename PixelIt>
inline bool PixelMerger<PixelIt>::Node::operator!=(const Node &other) const noexcept {
  return !(*this == other);
}

template <typename PixelIt>
inline PixelMerger<PixelIt>::PixelMerger(std::vector<PixelIt> heads, std::vector<PixelIt> tails) {
  assert(heads.size() == tails.size());
  for (std::size_t i = 0; i < heads.size(); ++i) {
    auto &first = heads[i];
    auto &last = tails[i];
    if (first != last) {
      _heads.emplace_back(std::move(first));
      _tails.emplace_back(std::move(last));
      _pqueue.emplace(Node{std::move(*_heads.back()++), _pqueue.size()});
    }
  }
}

template <typename PixelIt>
inline void PixelMerger<PixelIt>::replace_top_node(std::size_t i) {
  assert(this->_pqueue.top().i == i);
  this->_pqueue.pop();
  if (auto &it = this->_heads[i]; it != this->_tails[i]) {
    this->_pqueue.emplace(Node{*it++, i});
  }
}

template <typename PixelIt>
inline auto PixelMerger<PixelIt>::next() -> ThinPixel<N> {
  if (this->_pqueue.empty()) {
    return {};
  }

  auto current_node = this->_pqueue.top();
  this->replace_top_node(current_node.i);

  while (!this->_pqueue.empty()) {
    const auto next_node = this->_pqueue.top();
    if (next_node != current_node) {
      break;
    }
    current_node.pixel.count += next_node.pixel.count;
    this->replace_top_node(next_node.i);
  }
  return current_node.pixel;
}
}  // namespace internal

}  // namespace hictk
