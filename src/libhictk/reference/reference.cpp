// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/reference.hpp"

// clang-format off
#include "hictk/suppress_warnings.hpp"
HICTK_DISABLE_WARNING_PUSH
HICTK_DISABLE_WARNING_DEPRECATED_DECLARATIONS
#include <parallel_hashmap/phmap.h>
HICTK_DISABLE_WARNING_POP
// clang-format on

#include <fmt/format.h>
#include <fmt/std.h>

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
#include "hictk/numeric_utils.hpp"

namespace hictk {

Reference::Reference(std::initializer_list<Chromosome> chromosomes)
    : Reference(chromosomes.begin(), chromosomes.end()) {}

Reference Reference::from_chrom_sizes(const std::filesystem::path& path_to_chrom_sizes) {
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

auto Reference::begin() const -> const_iterator { return cbegin(); }
auto Reference::end() const -> const_iterator { return cend(); }
auto Reference::cbegin() const -> const_iterator { return _buff.cbegin(); }
auto Reference::cend() const -> const_iterator { return _buff.cend(); }

auto Reference::rbegin() const -> const_reverse_iterator { return rcbegin(); }
auto Reference::rend() const -> const_reverse_iterator { return rcend(); }
auto Reference::rcbegin() const -> const_reverse_iterator { return _buff.rbegin(); }
auto Reference::rcend() const -> const_reverse_iterator { return _buff.rend(); }

bool Reference::empty() const noexcept { return size() == 0; }
std::size_t Reference::size() const noexcept { return _buff.size(); }

auto Reference::find(std::uint32_t id) const -> const_iterator {
  if (static_cast<std::size_t>(id) > size()) {
    return end();
  }
  return _buff.begin() + static_cast<std::ptrdiff_t>(id);
}

auto Reference::find(std::string_view chrom_name) const -> const_iterator {
  auto it = _map.find(chrom_name);
  if (it == _map.end()) {
    return end();
  }

  return _buff.begin() + static_cast<std::ptrdiff_t>(it->second);
}

auto Reference::find(const Chromosome& chrom) const -> const_iterator {
  auto match = find(chrom.id());
  if (match != end() && *match != chrom) {
    match = end();
  }
  return match;
}

const Chromosome& Reference::at(std::uint32_t id) const {
  validate_chrom_id(id);
  return *find(id);
}

const Chromosome& Reference::at(std::string_view chrom_name) const {
  if (const auto match = find(chrom_name); match != end()) {
    return *match;
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

const Chromosome& Reference::operator[](std::uint32_t id) const noexcept {
  auto it = find(id);
  assert(it != end());
  return *it;
}
const Chromosome& Reference::operator[](std::string_view chrom_name) const noexcept {
  auto it = find(chrom_name);
  assert(it != end());
  return *it;
}

// NOLINTBEGIN(*-container-contains)
bool Reference::contains(std::uint32_t id) const { return find(id) != end(); }
bool Reference::contains(const Chromosome& chrom) const { return find(chrom) != end(); }
bool Reference::contains(std::string_view chrom_name) const { return find(chrom_name) != end(); }
// NOLINTEND(*-container-contains)

std::uint32_t Reference::get_id(std::string_view chrom_name) const {
  if (const auto match = find(chrom_name); match != end()) {
    return static_cast<std::uint32_t>(std::distance(begin(), match));
  }
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

bool Reference::operator==(const Reference& other) const {
  if (size() != other.size()) {
    return false;
  }
  return std::equal(_buff.begin(), _buff.end(), other.begin(),
                    [](const Chromosome& chrom1, const Chromosome& chrom2) {
                      // clang-format off
                      return chrom1.id() == chrom2.id()     &&
                             chrom1.name() == chrom2.name() &&
                             chrom1.size() == chrom2.size();
                      // clang-format on
                    });
}

bool Reference::operator!=(const Reference& other) const { return !(*this == other); }

const Chromosome& Reference::longest_chromosome() const {
  if (empty()) {
    throw std::runtime_error("longest_chromosome() was called on an empty Reference");
  }
  assert(_longest_chrom < _buff.size());
  return _buff[_longest_chrom];
}

const Chromosome& Reference::chromosome_with_longest_name() const {
  if (empty()) {
    throw std::runtime_error("chromosome_with_longest_name() was called on an empty Reference");
  }
  assert(_chrom_with_longest_name < _buff.size());
  return _buff[_chrom_with_longest_name];
}

Reference Reference::remove_ALL() const {
  std::vector<Chromosome> chroms{};
  std::copy_if(begin(), end(), std::back_inserter(chroms),
               [](const Chromosome& chrom) { return !chrom.is_all(); });

  return {chroms.begin(), chroms.end()};
}

Reference Reference::add_ALL(std::uint32_t scaling_factor) const {
  std::uint32_t all_size = 0;
  for (const auto& chrom : *this) {
    all_size += chrom.size() / scaling_factor;
  }

  std::vector chroms{Chromosome{0, "All", std::max(std::uint32_t{1}, all_size)}};
  std::copy_if(begin(), end(), std::back_inserter(chroms),
               [](const Chromosome& chrom) { return !chrom.is_all(); });

  return {chroms.begin(), chroms.end()};
}

void Reference::validate_chrom_id(std::uint32_t chrom_id) const {
  if (static_cast<std::size_t>(chrom_id) >= size()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }
}

auto Reference::construct_chrom_map(const ChromBuff& chroms) -> ChromMap {
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

std::size_t Reference::find_longest_chromosome(const ChromBuff& chroms) noexcept {
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

std::size_t Reference::find_chromosome_with_longest_name(const ChromBuff& chroms) noexcept {
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

std::vector<std::uint64_t> Reference::compute_size_prefix_sum(const ChromBuff& chroms) noexcept {
  std::vector<std::uint64_t> buff(chroms.size() + 2, 0);
  for (std::size_t i = 1; i < chroms.size() + 1; ++i) {
    buff[i] = buff[i - 1] + chroms[i - 1].size();
  }

  buff.back() = buff[chroms.size()] + 1;

  return buff;
}

void Reference::validate() const {
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
