// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <memory>
#include <string_view>
#include <vector>

namespace hictk {
class GenomicInterval;
class BinTable;
class Chromosome;
class Reference;

namespace cooler {

class Index {
  using ChromID = std::uint32_t;
  using OffsetVect = std::vector<std::uint64_t>;

  std::shared_ptr<const BinTable> _bins{};
  std::vector<OffsetVect> _idx{};
  std::size_t _size{};
  std::uint64_t _nnz{};

  static constexpr auto offset_not_set_value = std::numeric_limits<std::uint64_t>::max();

  // We pretend this is a map-like structure whose keys are used as index in the _idx vector
  using MapT = decltype(_idx);

 public:
  using key_type = ChromID;
  using value_type = typename MapT::value_type;
  using mapped_type = OffsetVect;
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
  explicit Index(std::shared_ptr<const BinTable> bins, std::uint64_t nnz = 0);

  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins() const noexcept;
  [[nodiscard]] std::shared_ptr<const BinTable> bins_ptr() const noexcept;

  [[nodiscard]] std::size_t num_chromosomes() const noexcept;
  [[nodiscard]] constexpr std::size_t size() const noexcept { return this->_size; }
  [[nodiscard]] std::size_t size(std::string_view chrom_name) const;
  [[nodiscard]] std::size_t size(std::uint32_t chrom_id) const;

  [[nodiscard]] constexpr bool empty() const noexcept { return this->size() == 0; }

  [[nodiscard]] std::uint32_t bin_size() const noexcept;

  [[nodiscard]] auto begin() const noexcept -> const_iterator;
  [[nodiscard]] auto end() const noexcept -> const_iterator;

  [[nodiscard]] auto cbegin() const noexcept -> const_iterator;
  [[nodiscard]] auto cend() const noexcept -> const_iterator;

  [[nodiscard]] auto at(std::string_view chrom_name) const -> const mapped_type&;
  [[nodiscard]] auto at(std::uint32_t chrom_id) const -> const mapped_type&;

  [[nodiscard]] std::uint64_t get_offset_by_bin_id(std::uint64_t bin_id) const;
  [[nodiscard]] std::uint64_t get_offset_by_pos(const Chromosome& chrom, std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t get_offset_by_pos(std::string_view chrom_name,
                                                std::uint32_t pos) const;
  [[nodiscard]] std::uint64_t get_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos) const;

  [[nodiscard]] std::uint64_t get_offset_by_row_idx(std::uint32_t chrom_id,
                                                    std::size_t row_idx) const;

  void set_offset_by_bin_id(std::uint64_t bin_id, std::uint64_t offset);

  void set_offset_by_pos(const Chromosome& chrom, std::uint32_t pos, std::uint64_t offset);
  void set_offset_by_pos(std::string_view chrom_name, std::uint32_t pos, std::uint64_t offset);
  void set_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos, std::uint64_t offset);

  void set_offset_by_row_idx(std::uint32_t chrom_id, std::size_t row_idx, std::uint64_t offset);

  void validate() const;

  [[nodiscard]] constexpr std::uint64_t nnz() const noexcept { return this->_nnz; }
  [[nodiscard]] std::uint64_t& nnz() noexcept;

  void compute_chrom_offsets(std::vector<std::uint64_t>& buff) const noexcept;
  [[nodiscard]] std::vector<std::uint64_t> compute_chrom_offsets() const;
  // Return first bin1_id corresponding to chrom_name/id
  [[nodiscard]] std::uint64_t chrom_to_bin1_offset(std::string_view chrom_name) const;
  [[nodiscard]] std::uint64_t chrom_to_bin1_offset(std::uint32_t chrom_id) const;

  void finalize(std::uint64_t nnz);

 private:
  [[nodiscard]] auto at(std::string_view chrom_name) -> mapped_type&;
  [[nodiscard]] auto at(std::uint32_t chrom_id) -> mapped_type&;

  void validate_chrom_id(std::uint32_t chrom_id) const;
  void validate(const Chromosome& chrom) const;

  [[nodiscard]] static auto init(const Reference& chroms, std::uint32_t bin_size) -> MapT;

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

#include "../../../index_impl.hpp"
