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
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/numeric_utils.hpp"

namespace hictk {

Reference::Reference(std::initializer_list<Chromosome> chromosomes)
    : Reference(chromosomes.begin(), chromosomes.end()) {}

[[nodiscard]] static std::vector<Chromosome> sanitize_chromosomes(
    std::vector<Chromosome> chromosomes, bool validate) {
  if (chromosomes.empty()) {
    return chromosomes;
  }

  if (!validate) {
    chromosomes.shrink_to_fit();
    return chromosomes;
  }

  if (!std::is_sorted(chromosomes.begin(), chromosomes.end())) {
    throw std::runtime_error("chromosomes are not sorted by ID");
  }

  phmap::flat_hash_set<std::uint32_t> ids(chromosomes.size());
  for (const auto& chrom : chromosomes) {
    if (chrom.size() == 0) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("chromosome {} has a size of 0"), chrom.name()));
    }
    ids.emplace(chrom.id());
  }

  if (ids.size() != chromosomes.size()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("found two or more chromosomes with the same ID")));
  }

  chromosomes.shrink_to_fit();
  return chromosomes;
}

Reference::Reference(ChromBuff chromosomes, bool validate)
    : _buff(sanitize_chromosomes(std::move(chromosomes), validate)),
      _map(construct_chrom_map(_buff)),
      _size_prefix_sum(compute_size_prefix_sum(_buff)),
      _longest_chrom(find_longest_chromosome(_buff)),
      _chrom_with_longest_name(find_chromosome_with_longest_name(_buff)) {
  if (empty()) {
    assert(_longest_chrom == 0);
    assert(_chrom_with_longest_name == 0);
    return;
  }

  assert(_longest_chrom < _buff.size());
  assert(_chrom_with_longest_name < _buff.size());
}

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
  const auto it = _map.find(chrom_name);
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

[[noreturn]] static void raise_chromosome_not_found(std::string_view chrom_name) {
  throw std::out_of_range(fmt::format(FMT_STRING("chromosome \"{}\" not found"), chrom_name));
}

const Chromosome& Reference::at(std::uint32_t id) const {
  validate_chrom_id(id);
  return *find(id);
}

const Chromosome& Reference::at(std::string_view chrom_name) const {
  if (const auto match = find(chrom_name); match != end()) {
    return *match;
  }
  raise_chromosome_not_found(chrom_name);
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
  raise_chromosome_not_found(chrom_name);
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
  if (HICTK_UNLIKELY(empty())) {
    throw std::runtime_error("longest_chromosome() was called on an empty Reference");
  }
  assert(_longest_chrom < _buff.size());
  return _buff[_longest_chrom];
}

const Chromosome& Reference::chromosome_with_longest_name() const {
  if (HICTK_UNLIKELY(empty())) {
    throw std::runtime_error("chromosome_with_longest_name() was called on an empty Reference");
  }
  assert(_chrom_with_longest_name < _buff.size());
  return _buff[_chrom_with_longest_name];
}

[[nodiscard]] static bool contains_ALL(const Reference& chroms) {
  return std::any_of(chroms.begin(), chroms.end(),
                     [](const Chromosome& chrom) { return chrom.is_all(); });
}

static void chrom_copy_helper(const Reference& chroms, std::vector<Chromosome>& buff) {
  for (const auto& chrom : chroms) {
    if (chrom.is_all()) {
      continue;
    }

    const auto chrom_id = static_cast<std::uint32_t>(buff.size());
    buff.emplace_back(chrom_id, chrom.name_ptr(), chrom.size());
  }
}

Reference Reference::remove_ALL() const {
  if (!contains_ALL(*this)) {
    return *this;
  }

  std::vector<Chromosome> chroms{};
  chroms.reserve(size() - 1);
  chrom_copy_helper(*this, chroms);

  return Reference{std::move(chroms), false};
}

[[nodiscard]] static std::uint32_t compute_ALL_size(const Reference& chroms,
                                                    std::uint32_t scaling_factor) noexcept {
  assert(scaling_factor != 0);
  const auto size =
      std::accumulate(chroms.begin(), chroms.end(), std::uint32_t{},
                      [scaling_factor](std::uint32_t accumulator, const Chromosome& chrom) {
                        return accumulator + (chrom.size() / scaling_factor);
                      });

  return std::max(std::uint32_t{1}, size);
}

Reference Reference::add_ALL(std::uint32_t scaling_factor) const {
  const auto all_size = compute_ALL_size(*this, scaling_factor);
  std::vector<Chromosome> chroms{};

  if (contains_ALL(*this)) {
    chroms.reserve(size());
  } else {
    chroms.reserve(size() + 1);
  }

  chroms.emplace_back(0, "All", all_size);
  chrom_copy_helper(*this, chroms);

  return Reference{std::move(chroms), false};
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

struct FindLongestChrom {};
struct FindLongestChromName {};

[[nodiscard]] static std::size_t get_size(const Chromosome& chrom, const FindLongestChrom&) {
  return chrom.size();
}

[[nodiscard]] static std::size_t get_size(const Chromosome& chrom, const FindLongestChromName&) {
  return chrom.name().size();
}

template <typename Strategy>
[[nodiscard]] static std::size_t find_longest_chromosome_impl(
    const std::vector<Chromosome>& chroms) noexcept {
  if (chroms.empty()) {
    return Chromosome{}.id();
  }

  std::size_t chrom_id = 0;
  for (std::size_t i = 1; i < chroms.size(); ++i) {
    const auto& chrom = chroms[i];
    if (chrom.is_all()) {
      continue;
    }

    const auto new_size = get_size(chrom, Strategy{});
    const auto max_size = get_size(chroms[chrom_id], Strategy{});

    if (new_size > max_size) {
      chrom_id = i;
    }
  }

  return chrom_id;
}

std::size_t Reference::find_longest_chromosome(const ChromBuff& chroms) noexcept {
  return find_longest_chromosome_impl<FindLongestChrom>(chroms);
}

std::size_t Reference::find_chromosome_with_longest_name(const ChromBuff& chroms) noexcept {
  return find_longest_chromosome_impl<FindLongestChromName>(chroms);
}

std::vector<std::uint64_t> Reference::compute_size_prefix_sum(const ChromBuff& chroms) noexcept {
  std::vector<std::uint64_t> buff(chroms.size() + 2, 0);
  for (std::size_t i = 1; i < chroms.size() + 1; ++i) {
    buff[i] = buff[i - 1] + chroms[i - 1].size();
  }

  buff.back() = buff[chroms.size()] + 1;

  return buff;
}

}  // namespace hictk
