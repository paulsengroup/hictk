// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <parallel_hashmap/phmap.h>

#include <cstdint>
#include <filesystem>
#include <initializer_list>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "hictk/chromosome.hpp"

namespace hictk {

class Reference {
  using ChromBuff = std::vector<Chromosome>;
  using ChromMap = phmap::flat_hash_map<std::string_view, std::size_t>;

  ChromBuff _buff{};
  ChromMap _map{};

  std::size_t _longest_chrom{Chromosome{}.id()};
  std::size_t _chrom_with_longest_name{Chromosome{}.id()};

 public:
  using value_type = typename ChromBuff::value_type;
  using size_type = typename ChromBuff::size_type;
  using difference_type = typename ChromBuff::difference_type;
  using allocator_type = typename ChromBuff::allocator_type;
  using reference = typename ChromBuff::const_reference;
  using const_reference = typename ChromBuff::const_reference;
  using pointer = typename ChromBuff::const_pointer;
  using const_pointer = typename ChromBuff::const_pointer;
  using iterator = typename ChromBuff::const_iterator;
  using const_iterator = typename ChromBuff::const_iterator;
  using reverse_iterator = typename ChromBuff::reverse_iterator;
  using const_reverse_iterator = typename ChromBuff::const_reverse_iterator;

  Reference() = default;

  template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
  Reference(ChromosomeNameIt first_chrom_name, ChromosomeNameIt last_chrom_name,
            ChromosomeSizeIt first_chrom_size);

  // Note: chromosome IDs are not preserved
  template <typename ChromosomeIt>
  Reference(ChromosomeIt first_chrom, ChromosomeIt last_chrom);
  Reference(std::initializer_list<Chromosome> chromosomes);
  [[nodiscard]] static Reference from_chrom_sizes(const std::filesystem::path& path_to_chrom_sizes);

  [[nodiscard]] auto begin() const -> const_iterator;
  [[nodiscard]] auto end() const -> const_iterator;
  [[nodiscard]] auto cbegin() const -> const_iterator;
  [[nodiscard]] auto cend() const -> const_iterator;

  [[nodiscard]] auto rbegin() const -> const_reverse_iterator;
  [[nodiscard]] auto rend() const -> const_reverse_iterator;
  [[nodiscard]] auto rcbegin() const -> const_reverse_iterator;
  [[nodiscard]] auto rcend() const -> const_reverse_iterator;

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] std::size_t size() const noexcept;

  [[nodiscard]] auto find(std::uint32_t id) const -> const_iterator;
  [[nodiscard]] auto find(std::string_view chrom_name) const -> const_iterator;
  [[nodiscard]] auto find(const Chromosome& chrom) const -> const_iterator;

  [[nodiscard]] const Chromosome& at(std::uint32_t id) const;
  [[nodiscard]] const Chromosome& at(std::string_view chrom_name) const;

  [[nodiscard]] const Chromosome& operator[](std::uint32_t id) const noexcept;
  [[nodiscard]] const Chromosome& operator[](std::string_view chrom_name) const noexcept;

  [[nodiscard]] bool contains(std::uint32_t id) const;
  [[nodiscard]] bool contains(const Chromosome& chrom) const;
  [[nodiscard]] bool contains(std::string_view chrom_name) const;

  [[nodiscard]] std::uint32_t get_id(std::string_view chrom_name) const;

  [[nodiscard]] bool operator==(const Reference& other) const;
  [[nodiscard]] bool operator!=(const Reference& other) const;

  // In case of ties, the first match is returned
  [[nodiscard]] const Chromosome& longest_chromosome() const;
  [[nodiscard]] const Chromosome& chromosome_with_longest_name() const;

 private:
  void validate_chrom_id(std::uint32_t chrom_id) const;

  template <typename ChromosomeNameIt, typename ChromosomeSizeIt>
  [[nodiscard]] static auto construct_chrom_buffer(ChromosomeNameIt first_chrom_name,
                                                   ChromosomeNameIt last_chrom_name,
                                                   ChromosomeSizeIt first_chrom_size) -> ChromBuff;

  [[nodiscard]] static auto construct_chrom_map(const ChromBuff& chroms) -> ChromMap;

  [[nodiscard]] static std::size_t find_longest_chromosome(const ChromBuff& chroms) noexcept;
  [[nodiscard]] static std::size_t find_chromosome_with_longest_name(
      const ChromBuff& chroms) noexcept;

  void validate() const;
};
}  // namespace hictk

#include "../../reference_impl.hpp"
