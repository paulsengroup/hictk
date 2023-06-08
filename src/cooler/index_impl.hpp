// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/fmt.hpp"

namespace hictk {

inline Index::Index(std::shared_ptr<const BinTable> bins, std::uint64_t nnz)
    : _bins(std::move(bins)),
      _idx(Index::init(_bins->chromosomes(), _bins->bin_size())),
      _nnz(nnz) {
  assert(this->bin_size() != 0);
  _size = std::accumulate(_idx.begin(), _idx.end(), std::size_t(0),
                          [&](std::size_t sum, const auto &it) { return sum + it.size(); });
}

inline const Reference &Index::chromosomes() const noexcept {
  assert(this->_bins);
  return this->_bins->chromosomes();
}

inline const BinTable &Index::bins() const noexcept {
  assert(this->_bins);
  return *this->_bins;
}

inline std::shared_ptr<const BinTable> Index::bins_ptr() const noexcept { return this->_bins; }

inline std::size_t Index::num_chromosomes() const noexcept {
  assert(this->_idx.size() == this->_bins->num_chromosomes());
  return this->_idx.size();
}

inline std::size_t Index::size(std::string_view chrom_name) const {
  const auto chrom_id = this->chromosomes().get_id(chrom_name);
  return this->size(chrom_id);
}

inline std::size_t Index::size(std::uint32_t chrom_id) const {
  this->validate_chrom_id(chrom_id);
  return this->at(chrom_id).size();
}

inline std::uint32_t Index::bin_size() const noexcept {
  assert(this->_bins);
  return this->_bins->bin_size();
}

inline auto Index::begin() const noexcept -> const_iterator { return iterator{this}; }
inline auto Index::end() const noexcept -> const_iterator {
  return iterator::make_end_iterator(this);
}

inline auto Index::cbegin() const noexcept -> const_iterator { return this->begin(); }
inline auto Index::cend() const noexcept -> const_iterator { return this->end(); }

inline auto Index::at(std::string_view chrom_name) const -> const mapped_type & {
  const auto chrom_id = this->chromosomes().get_id(chrom_name);
  return this->_idx.at(chrom_id);
}

inline auto Index::at(std::uint32_t chrom_id) -> mapped_type & {
  this->validate_chrom_id(chrom_id);
  return this->_idx.at(chrom_id);
}

inline auto Index::at(std::string_view chrom_name) -> mapped_type & {
  const auto chrom_id = this->chromosomes().get_id(chrom_name);
  return this->_idx.at(chrom_id);
}

inline auto Index::at(std::uint32_t chrom_id) const -> const mapped_type & {
  this->validate_chrom_id(chrom_id);
  return this->_idx.at(chrom_id);
}

inline std::uint64_t Index::get_offset_by_bin_id(std::uint64_t bin_id) const {
  if (bin_id == this->size()) {
    return this->_idx.back().back();
  }
  const auto &coords = this->_bins->at(bin_id);
  return this->get_offset_by_pos(coords.chrom(), coords.start());
}

inline std::uint64_t Index::get_offset_by_pos(const Chromosome &chrom, std::uint32_t pos) const {
  return this->get_offset_by_pos(chrom.name(), pos);
}

inline std::uint64_t Index::get_offset_by_pos(std::string_view chrom_name,
                                              std::uint32_t pos) const {
  const auto row_idx = pos / this->bin_size();
  return this->get_offset_by_row_idx(this->chromosomes().get_id(chrom_name), row_idx);
}

inline std::uint64_t Index::get_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos) const {
  const auto row_idx = pos / this->bin_size();
  return this->get_offset_by_row_idx(chrom_id, row_idx);
}

inline std::uint64_t Index::get_offset_by_row_idx(std::uint32_t chrom_id,
                                                  std::size_t row_idx) const {
  const auto &offsets = this->at(chrom_id);
  if (row_idx >= offsets.size()) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("invalid row_index {}: row maps outside of chromosome {}"), row_idx,
                    this->chromosomes().at(chrom_id)));
  }
  return offsets[row_idx];
}

inline void Index::set_offset_by_bin_id(std::uint64_t bin_id, std::uint64_t offset) {
  const auto &bin = this->_bins->at(bin_id);
  this->set_offset_by_pos(bin.chrom(), bin.start(), offset);
}

inline void Index::set_offset_by_pos(const Chromosome &chrom, std::uint32_t pos,
                                     std::uint64_t offset) {
  this->set_offset_by_pos(chrom.id(), pos, offset);
}

inline void Index::set_offset_by_pos(std::string_view chrom_name, std::uint32_t pos,
                                     std::uint64_t offset) {
  const auto row_idx = pos / this->bin_size();
  this->set_offset_by_row_idx(this->chromosomes().get_id(chrom_name), row_idx, offset);
}

inline void Index::set_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos,
                                     std::uint64_t offset) {
  const auto row_idx = pos / this->bin_size();
  this->set_offset_by_row_idx(chrom_id, row_idx, offset);
}

inline void Index::set_offset_by_row_idx(std::uint32_t chrom_id, std::size_t row_idx,
                                         std::uint64_t offset) {
  auto &offsets = this->at(chrom_id);
  assert(row_idx < offsets.size());
  offsets[row_idx] = offset;
}

inline void Index::validate() const {
  std::for_each(this->chromosomes().begin(), this->chromosomes().end(),
                [this](const Chromosome &chrom) { this->validate(chrom); });
}

inline std::uint64_t &Index::nnz() noexcept { return this->_nnz; }

inline std::vector<std::uint64_t> Index::compute_chrom_offsets() const {
  std::vector<std::uint64_t> buff(this->num_chromosomes());
  this->compute_chrom_offsets(buff);
  return buff;
}

inline std::uint64_t Index::chrom_to_bin1_offset(std::string_view chrom_name) const {
  return this->at(chrom_name).front();
}

inline std::uint64_t Index::chrom_to_bin1_offset(std::uint32_t chrom_id) const {
  return this->at(chrom_id).front();
}

inline void Index::finalize(std::uint64_t nnz) {
  this->_nnz = nnz;
  auto fill_value = nnz;

  std::for_each(this->_idx.rbegin(), this->_idx.rend(), [&](OffsetVect &offsets) {
    std::transform(offsets.rbegin(), offsets.rend(), offsets.rbegin(), [&fill_value](auto &offset) {
      if (offset == Index::offset_not_set_value) {
        return fill_value;
      }

      return fill_value = offset;
    });
  });
  assert(this->_idx[0][0] == 0 || this->_idx[0][0] == this->_idx[0][1]);
  this->_idx[0][0] = 0;
}

inline void Index::compute_chrom_offsets(std::vector<std::uint64_t> &buff) const noexcept {
  buff.resize(this->num_chromosomes() + 1);
  buff[0] = 0;

  std::transform(this->_idx.begin(), this->_idx.end(), buff.begin() + 1,
                 [offset = std::uint64_t(0)](const auto &it) mutable {
                   return offset += conditional_static_cast<std::uint64_t>(it.size());
                 });
}

inline void Index::validate_chrom_id(std::uint32_t chrom_id) const {
  if (static_cast<std::size_t>(chrom_id) >= this->num_chromosomes()) {
    throw std::out_of_range(fmt::format(FMT_STRING("chromosome with id {} not found"), chrom_id));
  }
}

inline auto Index::init(const Reference &chroms, std::uint32_t bin_size) -> MapT {
  assert(!chroms.empty());
  assert(bin_size != 0);

  MapT idx(chroms.size());
  std::transform(chroms.begin(), chroms.end(), idx.begin(), [&](const Chromosome &chrom) {
    const auto num_bins = (chrom.size() + bin_size - 1) / bin_size;
    return std::vector<std::uint64_t>(num_bins, Index::offset_not_set_value);
  });

  return idx;
}

inline void Index::validate(const Chromosome &chrom) const {
  try {
    const auto chrom_id = chrom.id();
    const auto &offsets = this->at(chrom_id);
    if (chrom_id == 0) {
      if (offsets.front() != 0) {
        throw std::runtime_error("first offset is not zero");
      }
    } else {
      const auto &prev_offsets = this->at(chrom_id - 1);
      if (offsets.front() < prev_offsets.back()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("offsets are not in ascending order: offset for "
                                   "bin {}:{}-{} should be >= {}, found {}"),
                        chrom.name(), 0, this->bin_size(), prev_offsets.back(), offsets.front()));
      }
    }

    if (const auto it = std::is_sorted_until(offsets.begin(), offsets.end()); it != offsets.end()) {
      const auto i = std::distance(offsets.begin(), it);
      assert(i != 0);
      throw std::runtime_error(
          fmt::format(FMT_STRING("offsets are not in ascending order: pixels/bin1_offset[{}]={} > "
                                 "pixels/bin1_offset[{}]={}\n"),
                      i - 1, *(it - 1), i, *(it)));
    }

    if (this->_nnz != 0) {
      auto match = std::find_if(offsets.begin(), offsets.end(),
                                [this](const auto offset) { return offset > this->_nnz; });
      if (match != offsets.end()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("invalid offset {}: offset is greater than nnz ({} > {})"),
                        *match, *match, this->_nnz));
      }
    }

  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("{} index is corrupted or incomplete: {}"), chrom.name(), e.what()));
  }
}

// NOLINTNEXTLINE
inline Index::iterator::iterator(const Index *idx) : _idx(idx), _chrom_id(0), _offset_idx(0) {
  assert(idx);
}

inline bool Index::iterator::operator==(const iterator &other) const noexcept {
  return this->_idx == other._idx && this->_chrom_id == other._chrom_id &&
         this->_offset_idx == other._offset_idx;
}

inline bool Index::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

inline auto Index::iterator::operator*() const -> value_type {
  assert(this->_idx);
  if (this->_chrom_id > this->last_chrom_id()) {
    return this->_idx->_nnz;
  }
  return this->get_offsets()[this->_offset_idx];
}

inline auto Index::iterator::operator++() -> iterator & {
  if (this->_chrom_id > this->last_chrom_id()) {
    return *this = make_end_iterator(this->_idx);
  }

  if (++this->_offset_idx >= this->get_offsets().size()) {
    if (++this->_chrom_id > this->last_chrom_id()) {
      return *this;  // Next dereference returns the index size
    }

    this->_offset_idx = 0;
  }

  return *this;
}

inline auto Index::iterator::operator++(int) -> iterator {
  auto it = *this;
  ++(*this);
  return it;
}

inline auto Index::iterator::make_end_iterator(const Index *idx) -> iterator {
  assert(idx);

  iterator it{};

  it._idx = idx;
  it._chrom_id = it.last_chrom_id() + 1;
  it._offset_idx = npos;

  return it;
}

inline std::uint32_t Index::iterator::last_chrom_id() const noexcept {
  assert(this->_idx);
  if (this->_idx->num_chromosomes() == 0) {
    return 0;
  }

  return static_cast<std::uint32_t>(this->_idx->num_chromosomes() - 1);
}

inline auto Index::iterator::get_offsets() const noexcept -> const OffsetVect & {
  assert(this->_chrom_id < static_cast<std::uint32_t>(this->_idx->size()));
  return this->_idx->_idx[static_cast<std::size_t>(this->_chrom_id)];
}

}  // namespace hictk
