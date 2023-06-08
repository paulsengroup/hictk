// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "coolerpp/common.hpp"
#include "coolerpp/internal/hash.hpp"

namespace coolerpp {

inline Chromosome::Chromosome(std::uint32_t id, std::string name_, std::uint32_t size_) noexcept
    : _name(std::move(name_)), _id(id), _size(size_) {
  assert(_id != (std::numeric_limits<std::uint32_t>::max)());
  assert(_size != 0);
}

constexpr Chromosome::operator bool() const noexcept { return this->id() != Chromosome::null_id; }

constexpr std::uint32_t Chromosome::id() const noexcept { return this->_id; }

inline std::string_view Chromosome::name() const noexcept { return this->_name; }

constexpr std::uint32_t Chromosome::size() const noexcept { return this->_size; }

constexpr bool Chromosome::operator<(const Chromosome& other) const noexcept {
  return this->id() < other.id();
}

constexpr bool Chromosome::operator>(const Chromosome& other) const noexcept {
  return this->id() > other.id();
}

constexpr bool Chromosome::operator<=(const Chromosome& other) const noexcept {
  return this->id() <= other.id();
}

constexpr bool Chromosome::operator>=(const Chromosome& other) const noexcept {
  return this->id() >= other.id();
}

inline bool Chromosome::operator==(const Chromosome& other) const noexcept {
  return this->id() == other.id() && this->name() == other.name() && this->size() == other.size();
}

inline bool Chromosome::operator!=(const Chromosome& other) const noexcept {
  return !(*this == other);
}

inline bool operator==(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name() == b_name;
}
inline bool operator!=(const Chromosome& a, std::string_view b_name) noexcept {
  return a.name() != b_name;
}

inline bool operator==(std::string_view a_name, const Chromosome& b) noexcept {
  return b == a_name;
}
inline bool operator!=(std::string_view a_name, const Chromosome& b) noexcept {
  return !(b == a_name);
}

constexpr bool operator<(const Chromosome& a, std::uint32_t b_id) noexcept { return a.id() < b_id; }
constexpr bool operator>(const Chromosome& a, std::uint32_t b_id) noexcept { return a.id() > b_id; }
constexpr bool operator<=(const Chromosome& a, std::uint32_t b_id) noexcept {
  return a.id() <= b_id;
}
constexpr bool operator>=(const Chromosome& a, std::uint32_t b_id) noexcept {
  return a.id() >= b_id;
}
constexpr bool operator==(const Chromosome& a, std::uint32_t b_id) noexcept {
  return a.id() == b_id;
}
constexpr bool operator!=(const Chromosome& a, std::uint32_t b_id) noexcept {
  return a.id() != b_id;
}

constexpr bool operator<(std::uint32_t a_id, const Chromosome& b) noexcept { return b > a_id; }
constexpr bool operator>(std::uint32_t a_id, const Chromosome& b) noexcept { return b < a_id; }
constexpr bool operator<=(std::uint32_t a_id, const Chromosome& b) noexcept { return b >= a_id; }
constexpr bool operator>=(std::uint32_t a_id, const Chromosome& b) noexcept { return b <= a_id; }
constexpr bool operator==(std::uint32_t a_id, const Chromosome& b) noexcept { return b == a_id; }
constexpr bool operator!=(std::uint32_t a_id, const Chromosome& b) noexcept { return b != a_id; }

template <typename ChromosomeIt>
inline ChromosomeSet::ChromosomeSet(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
    : _buff(first_chrom, last_chrom),
      _map(construct_chrom_map(_buff)),
      _longest_chrom(find_longest_chromosome(_buff)),
      _chrom_with_longest_name(find_chromosome_with_longest_name(_buff)) {
  this->validate();
}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
inline ChromosomeSet::ChromosomeSet(ChromosomeNameIt first_chrom_name,
                                    ChromosomeNameIt last_chrom_name,
                                    ChromosomeSizeIt first_chrom_size)
    : _buff(construct_chrom_buffer(first_chrom_name, last_chrom_name, first_chrom_size)),
      _map(construct_chrom_map(_buff)),
      _longest_chrom(find_longest_chromosome(_buff)),
      _chrom_with_longest_name(find_chromosome_with_longest_name(_buff)) {
  this->validate();
}

inline ChromosomeSet::ChromosomeSet(std::initializer_list<Chromosome> chromosomes)
    : ChromosomeSet(chromosomes.begin(), chromosomes.end()) {}

inline auto ChromosomeSet::begin() const -> const_iterator { return this->cbegin(); }
inline auto ChromosomeSet::end() const -> const_iterator { return this->cend(); }
inline auto ChromosomeSet::cbegin() const -> const_iterator { return this->_buff.cbegin(); }
inline auto ChromosomeSet::cend() const -> const_iterator { return this->_buff.cend(); }

inline auto ChromosomeSet::rbegin() const -> const_reverse_iterator { return this->rcbegin(); }
inline auto ChromosomeSet::rend() const -> const_reverse_iterator { return this->rcend(); }
inline auto ChromosomeSet::rcbegin() const -> const_reverse_iterator {
  return this->_buff.rbegin();
}
inline auto ChromosomeSet::rcend() const -> const_reverse_iterator { return this->_buff.rend(); }

inline bool ChromosomeSet::empty() const noexcept { return this->size() == 0; }
inline std::size_t ChromosomeSet::size() const noexcept { return this->_buff.size(); }

inline auto ChromosomeSet::find(std::uint32_t id) const -> const_iterator {
  if (static_cast<std::size_t>(id) > this->size()) {
    return this->end();
  }
  return this->_buff.begin() + static_cast<std::ptrdiff_t>(id);
}

inline auto ChromosomeSet::find(std::string_view chrom_name) const -> const_iterator {
  auto it = this->_map.find(chrom_name);
  if (it == this->_map.end()) {
    return this->end();
  }

  return this->_buff.begin() + static_cast<std::ptrdiff_t>(it->second);
}

inline auto ChromosomeSet::find(const Chromosome& chrom) const -> const_iterator {
  auto match = this->find(chrom.id());
  if (match != this->end() && *match != chrom) {
    match = this->end();
  }
  return match;
}

inline const Chromosome& ChromosomeSet::at(std::uint32_t id) const {
  this->validate_chrom_id(id);
  return *this->find(id);
}

inline const Chromosome& ChromosomeSet::at(std::string_view chrom_name) const {
  if (const auto match = this->find(chrom_name); match != this->end()) {
    return *match;
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

inline const Chromosome& ChromosomeSet::operator[](std::uint32_t id) const noexcept {
  auto it = this->find(id);
  assert(it != this->end());
  return *it;
}
inline const Chromosome& ChromosomeSet::operator[](std::string_view chrom_name) const noexcept {
  auto it = this->find(chrom_name);
  assert(it != this->end());
  return *it;
}

inline bool ChromosomeSet::contains(std::uint32_t id) const {
  return this->find(id) != this->end();
}
inline bool ChromosomeSet::contains(const Chromosome& chrom) const {
  return this->find(chrom) != this->end();
}
inline bool ChromosomeSet::contains(std::string_view chrom_name) const {
  return this->find(chrom_name) != this->end();
}

inline std::uint32_t ChromosomeSet::get_id(std::string_view chrom_name) const {
  if (const auto match = this->find(chrom_name); match != this->end()) {
    return static_cast<std::uint32_t>(std::distance(this->begin(), match));
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

inline bool ChromosomeSet::operator==(const ChromosomeSet& other) const {
  if (this->size() != other.size()) {
    return false;
  }
  return std::equal(this->_buff.begin(), this->_buff.end(), other.begin(),
                    [](const Chromosome& chrom1, const Chromosome& chrom2) {
                      return chrom1.id() == chrom2.id() && chrom1.name() == chrom2.name() &&
                             chrom1.size() == chrom2.size();
                    });
}

inline bool ChromosomeSet::operator!=(const ChromosomeSet& other) const {
  return !(*this == other);
}

inline const Chromosome& ChromosomeSet::longest_chromosome() const {
  if (this->empty()) {
    throw std::runtime_error("longest_chromosome() was called on an empty ChromosomeSet");
  }
  assert(this->_longest_chrom < this->_buff.size());
  return this->_buff[this->_longest_chrom];
}
inline const Chromosome& ChromosomeSet::chromosome_with_longest_name() const {
  if (this->empty()) {
    throw std::runtime_error("chromosome_with_longest_name() was called on an empty ChromosomeSet");
  }
  assert(this->_chrom_with_longest_name < this->_buff.size());
  return this->_buff[this->_chrom_with_longest_name];
}

inline void ChromosomeSet::validate_chrom_id(std::uint32_t chrom_id) const {
  if (static_cast<std::size_t>(chrom_id) >= this->size()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }
}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
inline auto ChromosomeSet::construct_chrom_buffer(ChromosomeNameIt first_chrom_name,
                                                  ChromosomeNameIt last_chrom_name,
                                                  ChromosomeSizeIt first_chrom_size) -> ChromBuff {
  ChromBuff buff{};
  while (first_chrom_name != last_chrom_name) {
    if (std::string_view{*first_chrom_name}.empty()) {
      throw std::runtime_error("found chromosome with empty name");
    }
    buff.emplace_back(static_cast<std::uint32_t>(buff.size()), std::string{*first_chrom_name},
                      conditional_static_cast<std::uint32_t>(*first_chrom_size));

    ++first_chrom_name;
    ++first_chrom_size;
  }
  return buff;
}

inline auto ChromosomeSet::construct_chrom_map(const ChromBuff& chroms) -> ChromMap {
  ChromMap buff(chroms.size());
  std::transform(chroms.begin(), chroms.end(), std::inserter(buff, buff.begin()),
                 [&](const auto& chrom) {
                   if (buff.contains(chrom.name())) {
                     throw std::runtime_error(fmt::format(
                         FMT_STRING("found multiple entries for chromosome \"{}\""), chrom.name()));
                   }
                   return std::make_pair(chrom.name(), static_cast<std::size_t>(chrom.id()));
                 });
  return buff;
}

inline std::size_t ChromosomeSet::find_longest_chromosome(const ChromBuff& chroms) noexcept {
  if (chroms.empty()) {
    return Chromosome{}.id();
  }

  const auto match = std::max_element(chroms.begin(), chroms.end(),
                                      [](const Chromosome& chrom1, const Chromosome& chrom2) {
                                        return chrom1.size() < chrom2.size();
                                      });

  return static_cast<std::size_t>(std::distance(chroms.begin(), match));
}
inline std::size_t ChromosomeSet::find_chromosome_with_longest_name(
    const ChromBuff& chroms) noexcept {
  if (chroms.empty()) {
    return Chromosome{}.id();
  }

  const auto match = std::max_element(chroms.begin(), chroms.end(),
                                      [](const Chromosome& chrom1, const Chromosome& chrom2) {
                                        return chrom1.name().size() < chrom2.name().size();
                                      });

  return static_cast<std::size_t>(std::distance(chroms.begin(), match));
}

inline void ChromosomeSet::validate() const {
  if (this->empty()) {
    return;
  }

  assert(this->_longest_chrom < this->_buff.size());
  assert(this->_chrom_with_longest_name < this->_buff.size());

  if (!std::is_sorted(this->_buff.begin(), this->_buff.end())) {
    throw std::runtime_error("chromosomes are not sorted by ID");
  }

  for (const auto& chrom : this->_buff) {
    if (chrom.size() == 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("chromosome {} has a size of 0"), chrom.name()));
    }
  }
}

}  // namespace coolerpp

inline std::size_t std::hash<coolerpp::Chromosome>::operator()(
    const coolerpp::Chromosome& c) const {
  return coolerpp::internal::hash_combine(0, c.id(), c.name(), c.size());
}

constexpr fmt::format_parse_context::iterator fmt::formatter<coolerpp::Chromosome>::parse(
    format_parse_context& ctx) {
  auto* it = ctx.begin();
  const auto* end = ctx.end();

  if (coolerpp::starts_with(ctx, "ucsc")) {
    this->presentation = Presentation::ucsc;
    it += std::string_view{"ucsc"}.size();
  } else if (coolerpp::starts_with(ctx, "tsv")) {
    this->presentation = Presentation::tsv;
    it += std::string_view{"tsv"}.size();
  }

  if (it != end && *it != '}') {
    throw fmt::format_error("invalid format");
  }

  return it;
}

inline fmt::format_context::iterator fmt::formatter<coolerpp::Chromosome>::format(
    const coolerpp::Chromosome& c, format_context& ctx) const {
  return this->presentation == Presentation::tsv
             ? fmt::format_to(ctx.out(), FMT_STRING("{}\t{}"), c.name(), c.size())
             : fmt::format_to(ctx.out(), FMT_STRING("{}:{}"), c.name(), c.size());
}
