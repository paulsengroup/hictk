// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/fmt/chromosome.hpp"

namespace hictk::cooler {

inline Index::Index(std::shared_ptr<const BinTable> bins,
                    const std::vector<std::uint64_t> &chrom_offsets, std::uint64_t nnz,
                    bool allocate)
    : _bins(std::move(bins)),
      _idx(Index::init(_bins->chromosomes(), *_bins, chrom_offsets, allocate)),
      _nnz(nnz) {}

inline const Reference &Index::chromosomes() const noexcept {
  assert(_bins);
  return _bins->chromosomes();
}

inline const BinTable &Index::bins() const noexcept {
  assert(_bins);
  return *_bins;
}

inline std::shared_ptr<const BinTable> Index::bins_ptr() const noexcept { return _bins; }

inline std::size_t Index::size() const noexcept { return bins().size(); }

inline std::size_t Index::size(std::string_view chrom_name) const {
  const auto chrom_id = chromosomes().get_id(chrom_name);
  return size(chrom_id);
}

inline std::size_t Index::size(std::uint32_t chrom_id) const { return at(chrom_id).size(); }

inline bool Index::empty() const noexcept { return size() == 0; }

inline bool Index::empty(std::uint32_t chrom_id) const noexcept {
  const auto &idx = at(chrom_id);
  return idx.front() == idx.back();
}

inline bool Index::empty(std::string_view chrom_name) const noexcept {
  return empty(chromosomes().at(chrom_name).id());
}

inline std::uint32_t Index::resolution() const noexcept {
  assert(_bins);
  return _bins->resolution();
}

inline auto Index::begin() const noexcept -> const_iterator { return iterator{this}; }
inline auto Index::end() const noexcept -> const_iterator {
  return iterator::make_end_iterator(this);
}

inline auto Index::cbegin() const noexcept -> const_iterator { return begin(); }
inline auto Index::cend() const noexcept -> const_iterator { return end(); }

inline auto Index::at(std::string_view chrom_name) const -> const mapped_type & {
  const auto chrom = chromosomes().at(chrom_name);
  return _idx.at(chrom);
}

inline auto Index::at(std::uint32_t chrom_id) -> mapped_type & { return _idx.at(chrom_id); }

inline auto Index::at(std::string_view chrom_name) -> mapped_type & {
  const auto chrom_id = chromosomes().get_id(chrom_name);
  return _idx.at(chrom_id);
}

inline auto Index::at(std::uint32_t chrom_id) const -> const mapped_type & {
  return _idx.at(chrom_id);
}

inline std::uint64_t Index::get_offset_by_bin_id(std::uint64_t bin_id) const {
  if (bin_id == size()) {
    return _nnz;
  }
  const auto &coords = _bins->at(bin_id);
  return get_offset_by_pos(coords.chrom(), coords.start());
}

inline std::uint64_t Index::get_offset_by_pos(const Chromosome &chrom, std::uint32_t pos) const {
  const auto row_idx = pos / resolution();
  return get_offset_by_row_idx(chrom.id(), row_idx);
}

inline std::uint64_t Index::get_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos) const {
  const auto row_idx = pos / resolution();
  return get_offset_by_row_idx(chrom_id, row_idx);
}

inline std::uint64_t Index::get_offset_by_row_idx(std::uint32_t chrom_id,
                                                  std::size_t row_idx) const {
  const auto &offsets = at(chrom_id);
  if (row_idx >= offsets.size()) {
    throw std::out_of_range(
        fmt::format(FMT_STRING("invalid row_index {}: row maps outside of chromosome {}"), row_idx,
                    chromosomes().at(chrom_id)));
  }
  return offsets[row_idx];
}

inline void Index::set(const Chromosome &chrom, OffsetVect offsets) {
  const auto [fist_bin, last_bin] = _bins->find_overlap(chrom, 0, chrom.size());
  const auto expected_size = static_cast<std::size_t>(std::distance(fist_bin, last_bin));
  if (offsets.size() != expected_size) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("expected index for {} to have size {}, found {}"), chrom,
                    expected_size, offsets.size()));
  }
  _idx.at(chrom) = std::move(offsets);
}

inline void Index::set_offset_by_bin(const Bin &bin, std::uint64_t offset) {
  set_offset_by_row_idx(bin.chrom().id(), conditional_static_cast<std::size_t>(bin.rel_id()),
                        offset);
}

inline void Index::set_offset_by_bin_id(std::uint64_t bin_id, std::uint64_t offset) {
  set_offset_by_bin(_bins->at(bin_id), offset);
}

inline void Index::set_offset_by_pos(const Chromosome &chrom, std::uint32_t pos,
                                     std::uint64_t offset) {
  set_offset_by_pos(chrom.id(), pos, offset);
}

inline void Index::set_offset_by_pos(std::uint32_t chrom_id, std::uint32_t pos,
                                     std::uint64_t offset) {
  set_offset_by_bin(_bins->at(chrom_id, pos), offset);
}

inline void Index::set_offset_by_row_idx(std::uint32_t chrom_id, std::size_t row_idx,
                                         std::uint64_t offset) {
  auto &offsets = at(chrom_id);
  assert(row_idx < offsets.size());
  offsets[row_idx] = offset;
}

inline void Index::validate() const {
  std::for_each(chromosomes().begin(), chromosomes().end(),
                [this](const Chromosome &chrom) { validate(chrom); });
}

constexpr std::uint64_t Index::nnz() const noexcept { return _nnz; }
inline void Index::set_nnz(std::uint64_t n) noexcept { _nnz = n; }

inline std::vector<std::uint64_t> Index::compute_chrom_offsets() const {
  std::vector<std::uint64_t> buff(chromosomes().size());
  compute_chrom_offsets(buff);
  return buff;
}

inline std::uint64_t Index::chrom_to_bin1_offset(std::string_view chrom_name) const {
  return at(chrom_name).front();
}

inline std::uint64_t Index::chrom_to_bin1_offset(std::uint32_t chrom_id) const {
  return at(chrom_id).front();
}

inline void Index::finalize(std::uint64_t nnz) {
  _nnz = nnz;
  auto fill_value = nnz;

  std::for_each(_idx.rbegin(), _idx.rend(), [&](auto &it) {
    auto &offsets = it.second;
    std::transform(offsets.rbegin(), offsets.rend(), offsets.rbegin(), [&fill_value](auto &offset) {
      if (offset == Index::offset_not_set_value) {
        return fill_value;
      }

      return fill_value = offset;
    });
  });
  _idx.begin()->second.front() = 0;
}

inline void Index::compute_chrom_offsets(std::vector<std::uint64_t> &buff) const noexcept {
  buff.resize(chromosomes().size() + 1);
  buff[0] = 0;

  std::transform(_idx.begin(), _idx.end(), buff.begin() + 1,
                 [offset = std::uint64_t(0)](const auto &it) mutable {
                   return offset += conditional_static_cast<std::uint64_t>(it.second.size());
                 });
}

inline auto Index::init(const Reference &chroms, const BinTable &bins,
                        const std::vector<std::uint64_t> &chrom_offsets, bool allocate) -> MapT {
  assert(!chroms.empty());
  assert(chrom_offsets.empty() || chroms.size() + 1 == chrom_offsets.size());
  MapT idx{};
  for (std::uint32_t i = 0; i < chroms.size(); ++i) {
    const auto &chrom = chroms.at(i);
    const auto [first_bin, last_bin] = bins.find_overlap(chrom, 0, chrom.size());
    const auto num_bins = static_cast<std::size_t>(std::distance(first_bin, last_bin));
    auto node =
        idx.emplace(chrom, OffsetVect(allocate ? num_bins : 1, Index::offset_not_set_value));

    if (!chrom_offsets.empty()) {
      const auto &offset = chrom_offsets.at(i);
      node.first->second.front() = offset;
    }
  }
  return idx;
}

inline void Index::validate(const Chromosome &chrom) const {
  try {
    const auto chrom_id = chrom.id();
    const auto &offsets = at(chrom_id);
    if (offsets.empty()) {
      throw std::runtime_error("offset vector is empty");
    }
    if (chrom_id == 0) {
      if (offsets.front() != 0) {
        throw std::runtime_error("first offset is not zero");
      }
    } else {
      const auto &prev_offsets = at(chrom_id - 1);
      if (offsets.front() < prev_offsets.back()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("offsets are not in ascending order: offset for "
                                   "bin {}:{}-{} should be >= {}, found {}"),
                        chrom.name(), 0, resolution(), prev_offsets.back(), offsets.front()));
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

    if (_nnz != 0) {
      auto match = std::find_if(offsets.begin(), offsets.end(),
                                [this](const auto offset) { return offset > _nnz; });
      if (match != offsets.end()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("invalid offset {}: offset is greater than nnz ({} > {})"),
                        *match, *match, _nnz));
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
  return _idx == other._idx && _chrom_id == other._chrom_id && _offset_idx == other._offset_idx;
}

inline bool Index::iterator::operator!=(const iterator &other) const noexcept {
  return !(*this == other);
}

inline auto Index::iterator::operator*() const -> value_type {
  assert(_idx);
  if (_chrom_id > last_chrom_id()) {
    return _idx->_nnz;
  }
  return get_offsets()[_offset_idx];
}

inline auto Index::iterator::operator++() -> iterator & {
  if (_chrom_id > last_chrom_id()) {
    return *this = make_end_iterator(_idx);
  }

  if (++_offset_idx >= get_offsets().size()) {
    if (++_chrom_id > last_chrom_id()) {
      return *this;  // Next dereference returns the index size
    }

    _offset_idx = 0;
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
  assert(_idx);
  if (_idx->size() == 0) {
    return 0;
  }

  return static_cast<std::uint32_t>(_idx->chromosomes().size() - 1);
}

inline auto Index::iterator::get_offsets() const noexcept -> const OffsetVect & {
  assert(_chrom_id < static_cast<std::uint32_t>(_idx->size()));
  return _idx->_idx.at(_chrom_id);
}

}  // namespace hictk::cooler
