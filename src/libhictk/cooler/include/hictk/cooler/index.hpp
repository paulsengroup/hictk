// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <parallel_hashmap/btree.h>

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <memory>
#include <string_view>
#include <vector>

#include "hictk/chromosome.hpp"

namespace hictk {
class GenomicInterval;
class BinTable;
class Reference;

namespace cooler {

class Index {
  using OffsetVect = std::vector<std::uint64_t>;
  using MapT = phmap::btree_map<Chromosome, OffsetVect, ChromosomeCmp>;

  std::shared_ptr<const BinTable> _bins{};
  MapT _idx{};
  std::uint64_t _nnz{};

  static constexpr auto offset_not_set_value = std::numeric_limits<std::uint64_t>::max();

 public:
  using key_type = typename MapT::key_type;
  using value_type = typename MapT::value_type;
  using mapped_type = typename MapT::mapped_type;
  using size_type = typename MapT::size_type;
  using difference_type = typename MapT::difference_type;
  using allocator_type = typename MapT::allocator_type;
  using reference = typename MapT::const_reference;
  using const_reference = typename MapT::const_reference;
  using pointer = typename MapT::const_pointer;
  using const_pointer = typename MapT::const_pointer;

  class iterator;
  using const_iterator = iterator;

  Index() = default;
  explicit Index(std::shared_ptr<const BinTable> bins,
                 const std::vector<std::uint64_t>& chrom_offsets = {}, std::uint64_t nnz = 0,
                 bool allocate = true);

  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;

  [[nodiscard]] std::size_t size() const noexcept;
  [[nodiscard]] std::size_t size(std::string_view chrom_name) const;
  [[nodiscard]] std::size_t size(std::uint32_t chrom_id) const;

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] bool empty(std::string_view chrom_name) const noexcept;
  [[nodiscard]] bool empty(std::uint32_t chrom_id) const noexcept;

  [[nodiscard]] std::uint32_t bin_size() const noexcept;

  [[nodiscard]] auto begin() const noexcept -> const_iterator;
  [[nodiscard]] auto end() const noexcept -> const_iterator;

  [[nodiscard]] auto cbegin() const noexcept -> const_iterator;
  [[nodiscard]] auto cend() const noexcept -> const_iterator;

  [[nodiscard]] auto at(std::string_view chrom_name) const -> const mapped_type&;
  [[nodiscard]] auto at(std::uint32_t chrom_id) const -> const mapped_type&;

  [[nodiscard]] std::uint64_t get_offset_by_bin_id(std::uint64_t bin_id) const;
  [[nodiscard]] std::uint64_t get_offset_by_pos(const Chromosome& chrom, std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t get_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos) const;

  [[nodiscard]] std::uint64_t get_offset_by_row_idx(std::uint32_t chrom_id,
                                                    std::size_t row_idx) const;

  void set(const Chromosome& chrom, OffsetVect offsets);
  void set_offset_by_bin_id(std::uint64_t bin_id, std::uint64_t offset);

  void set_offset_by_pos(const Chromosome& chrom, std::uint32_t pos, std::uint64_t offset);
  void set_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos, std::uint64_t offset);

  void set_offset_by_row_idx(std::uint32_t chrom_id, std::size_t row_idx, std::uint64_t offset);

  void validate() const;
  void validate(const Chromosome& chrom) const;

  [[nodiscard]] constexpr std::uint64_t nnz() const noexcept;
  void set_nnz(std::uint64_t n) noexcept;

  void compute_chrom_offsets(std::vector<std::uint64_t>& buff) const noexcept;
  [[nodiscard]] std::vector<std::uint64_t> compute_chrom_offsets() const;
  // Return first bin1_id corresponding to chrom_name/id
  [[nodiscard]] std::uint64_t chrom_to_bin1_offset(std::string_view chrom_name) const;
  [[nodiscard]] std::uint64_t chrom_to_bin1_offset(std::uint32_t chrom_id) const;

  void finalize(std::uint64_t nnz);

 private:
  [[nodiscard]] auto at(std::string_view chrom_name) -> mapped_type&;
  [[nodiscard]] auto at(std::uint32_t chrom_id) -> mapped_type&;

  [[nodiscard]] static auto init(const Reference& chroms, const BinTable& bins,
                                 const std::vector<std::uint64_t>& chrom_offsets, bool allocate)
      -> MapT;

 public:
  class iterator {
    friend Index;

    static constexpr auto npos = std::numeric_limits<std::size_t>::max();

    const Index* _idx{};
    std::uint32_t _chrom_id{};
    std::size_t _offset_idx{npos};

    iterator() = default;
    explicit iterator(const Index* idx);

   public:
    using difference_type = std::ptrdiff_t;
    using value_type = std::uint64_t;
    using pointer = value_type*;
    using reference = value_type&;
    using iterator_category = std::forward_iterator_tag;

    [[nodiscard]] bool operator==(const iterator& other) const noexcept;
    [[nodiscard]] bool operator!=(const iterator& other) const noexcept;
    [[nodiscard]] auto operator*() const -> value_type;
    auto operator++() -> iterator&;
    auto operator++(int) -> iterator;

   private:
    [[nodiscard]] std::uint32_t last_chrom_id() const noexcept;
    [[nodiscard]] auto get_offsets() const noexcept -> const OffsetVect&;
    [[nodiscard]] static auto make_end_iterator(const Index* idx) -> iterator;
  };
};

}  // namespace cooler

}  // namespace hictk

#include "./impl/index_impl.hpp"
