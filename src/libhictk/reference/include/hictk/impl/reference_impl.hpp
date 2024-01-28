// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <fmt/std.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/numeric_utils.hpp"

namespace hictk {

template <typename ChromosomeIt>
inline Reference::Reference(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
    : _buff(construct_chrom_buffer(first_chrom, last_chrom)),
      _map(construct_chrom_map(_buff)),
      _size_prefix_sum(compute_size_prefix_sum(_buff)),
      _longest_chrom(find_longest_chromosome(_buff)),
      _chrom_with_longest_name(find_chromosome_with_longest_name(_buff)) {
  validate();
}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
inline Reference::Reference(ChromosomeNameIt first_chrom_name, ChromosomeNameIt last_chrom_name,
                            ChromosomeSizeIt first_chrom_size)
    : _buff(construct_chrom_buffer(first_chrom_name, last_chrom_name, first_chrom_size)),
      _map(construct_chrom_map(_buff)),
      _size_prefix_sum(compute_size_prefix_sum(_buff)),
      _longest_chrom(find_longest_chromosome(_buff)),
      _chrom_with_longest_name(find_chromosome_with_longest_name(_buff)) {
  validate();
}

inline Reference::Reference(std::initializer_list<Chromosome> chromosomes)
    : Reference(chromosomes.begin(), chromosomes.end()) {}

inline Reference Reference::from_chrom_sizes(const std::filesystem::path& path_to_chrom_sizes) {
  try {
    std::string line;
    std::vector<std::string> chrom_names{};
    std::vector<std::uint32_t> chrom_sizes{};

    std::ifstream f(path_to_chrom_sizes);
    while (std::getline(f, line)) {
      const auto delim_pos = line.find('\t');

      chrom_names.emplace_back(line.substr(0, delim_pos));
      chrom_sizes.emplace_back(
          internal::parse_numeric_or_throw<std::uint32_t>(line.substr(delim_pos + 1)));
    }

    return {chrom_names.begin(), chrom_names.end(), chrom_sizes.begin()};
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while importing chromosomes from {}: {}"),
                    path_to_chrom_sizes, e.what()));
  }
}

inline auto Reference::begin() const -> const_iterator { return cbegin(); }
inline auto Reference::end() const -> const_iterator { return cend(); }
inline auto Reference::cbegin() const -> const_iterator { return _buff.cbegin(); }
inline auto Reference::cend() const -> const_iterator { return _buff.cend(); }

inline auto Reference::rbegin() const -> const_reverse_iterator { return rcbegin(); }
inline auto Reference::rend() const -> const_reverse_iterator { return rcend(); }
inline auto Reference::rcbegin() const -> const_reverse_iterator { return _buff.rbegin(); }
inline auto Reference::rcend() const -> const_reverse_iterator { return _buff.rend(); }

inline bool Reference::empty() const noexcept { return size() == 0; }
inline std::size_t Reference::size() const noexcept { return _buff.size(); }

inline auto Reference::find(std::uint32_t id) const -> const_iterator {
  if (static_cast<std::size_t>(id) > size()) {
    return end();
  }
  return _buff.begin() + static_cast<std::ptrdiff_t>(id);
}

inline auto Reference::find(std::string_view chrom_name) const -> const_iterator {
  auto it = _map.find(chrom_name);
  if (it == _map.end()) {
    return end();
  }

  return _buff.begin() + static_cast<std::ptrdiff_t>(it->second);
}

inline auto Reference::find(const Chromosome& chrom) const -> const_iterator {
  auto match = find(chrom.id());
  if (match != end() && *match != chrom) {
    match = end();
  }
  return match;
}

inline const Chromosome& Reference::at(std::uint32_t id) const {
  validate_chrom_id(id);
  return *find(id);
}

inline const Chromosome& Reference::at(std::string_view chrom_name) const {
  if (const auto match = find(chrom_name); match != end()) {
    return *match;
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

inline const Chromosome& Reference::operator[](std::uint32_t id) const noexcept {
  auto it = find(id);
  assert(it != end());
  return *it;
}
inline const Chromosome& Reference::operator[](std::string_view chrom_name) const noexcept {
  auto it = find(chrom_name);
  assert(it != end());
  return *it;
}

inline bool Reference::contains(std::uint32_t id) const { return find(id) != end(); }
inline bool Reference::contains(const Chromosome& chrom) const { return find(chrom) != end(); }
inline bool Reference::contains(std::string_view chrom_name) const {
  return find(chrom_name) != end();
}

inline std::uint32_t Reference::get_id(std::string_view chrom_name) const {
  if (const auto match = find(chrom_name); match != end()) {
    return static_cast<std::uint32_t>(std::distance(begin(), match));
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

inline bool Reference::operator==(const Reference& other) const {
  if (size() != other.size()) {
    return false;
  }
  return std::equal(_buff.begin(), _buff.end(), other.begin(),
                    [](const Chromosome& chrom1, const Chromosome& chrom2) {
                      return chrom1.id() == chrom2.id() && chrom1.name() == chrom2.name() &&
                             chrom1.size() == chrom2.size();
                    });
}

inline bool Reference::operator!=(const Reference& other) const { return !(*this == other); }

constexpr const std::vector<std::uint64_t>& Reference::chrom_size_prefix_sum() const noexcept {
  return _size_prefix_sum;
}

inline const Chromosome& Reference::longest_chromosome() const {
  if (empty()) {
    throw std::runtime_error("longest_chromosome() was called on an empty Reference");
  }
  assert(_longest_chrom < _buff.size());
  return _buff[_longest_chrom];
}
inline const Chromosome& Reference::chromosome_with_longest_name() const {
  if (empty()) {
    throw std::runtime_error("chromosome_with_longest_name() was called on an empty Reference");
  }
  assert(_chrom_with_longest_name < _buff.size());
  return _buff[_chrom_with_longest_name];
}

inline Reference Reference::remove_ALL() const {
  std::vector<Chromosome> chroms{};
  std::copy_if(begin(), end(), std::back_inserter(chroms),
               [](const Chromosome& chrom) { return !chrom.is_all(); });

  return {chroms.begin(), chroms.end()};
}

inline Reference Reference::add_ALL(std::uint32_t scaling_factor) const {
  std::uint32_t all_size = 0;
  for (const auto& chrom : *this) {
    all_size += chrom.size() / scaling_factor;
  }

  std::vector<Chromosome> chroms{Chromosome{0, "All", all_size}};
  std::copy_if(begin(), end(), std::back_inserter(chroms),
               [](const Chromosome& chrom) { return !chrom.is_all(); });

  return {chroms.begin(), chroms.end()};
}

inline void Reference::validate_chrom_id(std::uint32_t chrom_id) const {
  if (static_cast<std::size_t>(chrom_id) >= size()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }
}

template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
inline auto Reference::construct_chrom_buffer(ChromosomeNameIt first_chrom_name,
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

template <typename ChromosomeIt>
inline auto Reference::construct_chrom_buffer(ChromosomeIt first_chrom, ChromosomeIt last_chrom)
    -> ChromBuff {
  std::vector<std::string> chrom_names{};
  std::vector<std::uint32_t> chrom_sizes{};

  std::for_each(first_chrom, last_chrom, [&](const Chromosome& chrom) {
    chrom_names.emplace_back(std::string{chrom.name()});
    chrom_sizes.emplace_back(chrom.size());
  });

  return construct_chrom_buffer(chrom_names.begin(), chrom_names.end(), chrom_sizes.begin());
}

inline auto Reference::construct_chrom_map(const ChromBuff& chroms) -> ChromMap {
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

inline std::size_t Reference::find_longest_chromosome(const ChromBuff& chroms) noexcept {
  if (chroms.empty()) {
    return Chromosome{}.id();
  }

  std::uint32_t max_length = 0;
  std::size_t i = Chromosome{}.id();
  for (std::size_t j = 0; j < chroms.size(); ++j) {
    const auto& chrom = chroms[j];
    if (chrom.is_all()) {
      continue;
    }
    if (chrom.size() > max_length) {
      max_length = chrom.size();
      i = j;
    }
  }

  return i;
}
inline std::size_t Reference::find_chromosome_with_longest_name(const ChromBuff& chroms) noexcept {
  if (chroms.empty()) {
    return Chromosome{}.id();
  }

  std::size_t max_length = 0;
  std::size_t i = Chromosome{}.id();
  for (std::size_t j = 0; j < chroms.size(); ++j) {
    const auto& chrom = chroms[j];
    if (chrom.is_all()) {
      continue;
    }
    if (chrom.name().size() > max_length) {
      max_length = chrom.name().size();
      i = j;
    }
  }

  return i;
}

inline std::vector<std::uint64_t> Reference::compute_size_prefix_sum(
    const ChromBuff& chroms) noexcept {
  std::vector<std::uint64_t> buff(chroms.size() + 2, 0);
  for (std::size_t i = 1; i < chroms.size() + 1; ++i) {
    buff[i] = buff[i - 1] + chroms[i - 1].size();
  }

  buff.back() = buff[chroms.size()] + 1;

  return buff;
}

inline void Reference::validate() const {
  if (empty()) {
    return;
  }

  assert(_longest_chrom < _buff.size());
  assert(_chrom_with_longest_name < _buff.size());

  if (!std::is_sorted(_buff.begin(), _buff.end())) {
    throw std::runtime_error("chromosomes are not sorted by ID");
  }

  phmap::flat_hash_set<std::uint32_t> ids{};
  for (const auto& chrom : _buff) {
    if (chrom.size() == 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("chromosome {} has a size of 0"), chrom.name()));
    }
    ids.emplace(chrom.id());
  }

  if (ids.size() != _buff.size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("found two or more chromosomes with the same ID")));
  }
}

}  // namespace hictk
