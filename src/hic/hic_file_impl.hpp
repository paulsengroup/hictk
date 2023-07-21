// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "hictk/balancing/weights.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/footer.hpp"

namespace hictk::hic {

inline HiCFile::HiCFile(std::string url_, std::uint32_t resolution_, MatrixType type_,
                        MatrixUnit unit_, std::uint64_t block_cache_capacity)
    : _fs(std::make_shared<internal::HiCFileReader>(std::move(url_))),
      _type(type_),
      _unit(unit_),
      _block_cache(std::make_shared<internal::BlockCache>(block_cache_capacity)),
      _weight_cache(std::make_shared<internal::WeightCache>()),
      _bins(std::make_shared<const BinTable>(_fs->header().chromosomes, resolution_)) {
  if (!has_resolution(resolution())) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("file {} does not have interactions for resolution {}"), url(), resolution()));
  }

  if (block_cache_capacity == 0) {
    optimize_cache_size();
  }
}

inline HiCFile& HiCFile::open(std::string url_, std::uint32_t resolution_, MatrixType type_,
                              MatrixUnit unit_, std::uint64_t block_cache_capacity) {
  if (_fs->url() == url_ && resolution() == resolution_ && _type == type_ && _unit == unit_) {
    _block_cache->set_capacity(block_cache_capacity, false);
    return *this;
  }

  const auto prev_block_cache_capacity = _block_cache->capacity_bytes();
  *this = HiCFile(url_, resolution_, type_, unit_, block_cache_capacity);

  if (_block_cache->capacity_bytes() < prev_block_cache_capacity) {
    _block_cache->set_capacity(prev_block_cache_capacity);
  }
  return *this;
}

inline HiCFile& HiCFile::open(std::uint32_t resolution_, MatrixType type_, MatrixUnit unit_,
                              std::uint64_t block_cache_capacity) {
  return open(url(), resolution_, type_, unit_, block_cache_capacity);
}

inline bool HiCFile::has_resolution(std::uint32_t resolution) const {
  const auto match = std::find(avail_resolutions().begin(), avail_resolutions().end(), resolution);
  return match != avail_resolutions().end();
}

inline const std::string& HiCFile::url() const noexcept { return _fs->url(); }

inline const std::string& HiCFile::name() const noexcept { return url(); }

inline std::int32_t HiCFile::version() const noexcept { return _fs->version(); }

inline const BinTable& HiCFile::bins() const noexcept {
  assert(_bins);
  return *_bins;
}

inline std::shared_ptr<const BinTable> HiCFile::bins_ptr() const noexcept { return _bins; }

inline const Reference& HiCFile::chromosomes() const noexcept { return bins().chromosomes(); }

inline const std::string& HiCFile::assembly() const noexcept { return _fs->header().genomeID; }

inline const std::vector<std::uint32_t>& HiCFile::avail_resolutions() const noexcept {
  return _fs->header().resolutions;
}

inline std::uint32_t HiCFile::resolution() const noexcept { return _bins->bin_size(); }

inline std::shared_ptr<const internal::HiCFooter> HiCFile::get_footer(
    const Chromosome& chrom1, const Chromosome& chrom2, MatrixType matrix_type,
    NormalizationMethod norm, MatrixUnit unit, std::uint32_t resolution) const {
  const internal::HiCFooterMetadata metadata{url(),      matrix_type, norm,  unit,
                                             resolution, chrom1,      chrom2};
  auto it = _footers.find(metadata);
  if (it != _footers.end()) {
    return *it;
  }

  auto weights1 = _weight_cache->find_or_emplace(chrom1, norm);
  auto weights2 = _weight_cache->find_or_emplace(chrom2, norm);

  auto [node, _] = _footers.emplace(_fs->read_footer(chrom1.id(), chrom2.id(), matrix_type, norm,
                                                     unit, resolution, weights1, weights2));

  return *node;
}

inline PixelSelectorAll HiCFile::fetch(NormalizationMethod norm) const {
  std::vector<PixelSelector> selectors;

  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto& chrom1 = chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto& chrom2 = chromosomes().at(chrom2_id);
      if (chrom2.is_all()) {
        continue;
      }
      selectors.emplace_back(fetch(chrom1.name(), chrom2.name(), norm));
    }
  }

  return PixelSelectorAll{std::move(selectors)};
}

inline PixelSelector HiCFile::fetch(std::string_view query, NormalizationMethod norm,
                                    QUERY_TYPE query_type) const {
  const auto gi = query_type == QUERY_TYPE::BED
                      ? GenomicInterval::parse_bed(chromosomes(), query)
                      : GenomicInterval::parse_ucsc(chromosomes(), std::string{query});

  return fetch(gi.chrom(), gi.start(), gi.end(), gi.chrom(), gi.start(), gi.end(), norm);
}

inline PixelSelector HiCFile::fetch(std::string_view chrom_name, std::uint32_t start,
                                    std::uint32_t end, NormalizationMethod norm) const {
  return fetch(chrom_name, start, end, chrom_name, start, end, norm);
}

inline PixelSelector HiCFile::fetch(std::string_view range1, std::string_view range2,
                                    NormalizationMethod norm, QUERY_TYPE query_type) const {
  const auto gi1 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(chromosomes(), range1)
                       : GenomicInterval::parse_ucsc(chromosomes(), std::string{range1});

  const auto gi2 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(chromosomes(), range2)
                       : GenomicInterval::parse_ucsc(chromosomes(), std::string{range2});

  return fetch(gi1.chrom(), gi1.start(), gi1.end(), gi2.chrom(), gi2.start(), gi2.end(), norm);
}

inline PixelSelector HiCFile::fetch(std::string_view chrom1_name, std::uint32_t start1,
                                    std::uint32_t end1, std::string_view chrom2_name,
                                    std::uint32_t start2, std::uint32_t end2,
                                    NormalizationMethod norm) const {
  return fetch(chromosomes().at(chrom1_name), start1, end1, chromosomes().at(chrom2_name), start2,
               end2, norm);
}

inline PixelSelector HiCFile::fetch(const Chromosome& chrom1, std::uint32_t start1,
                                    std::uint32_t end1, const Chromosome& chrom2,
                                    std::uint32_t start2, std::uint32_t end2,
                                    NormalizationMethod norm) const {
  if (chrom1 > chrom2) {
    throw std::runtime_error(
        "Query overlaps the lower-triangle of the matrix. This is currently not supported.");
  }

  if (_type == MatrixType::expected && norm != NormalizationMethod::NONE) {
    throw std::logic_error(fmt::format(
        FMT_STRING("matrix type {} is incompatible with normalization method {}"), _type, norm));
  }

  const PixelCoordinates coord1 = {_bins->at(chrom1, start1), _bins->at(chrom1, end1 - 1)};
  const PixelCoordinates coord2 = {_bins->at(chrom2, start2), _bins->at(chrom2, end2 - 1)};

  return {_fs,          get_footer(chrom1, chrom2, _type, norm, _unit, resolution()),
          _block_cache, _bins,
          coord1,       coord2};
}

inline std::size_t HiCFile::num_cached_footers() const noexcept { return _footers.size(); }

inline void HiCFile::purge_footer_cache() { _footers.clear(); }

inline double HiCFile::block_cache_hit_rate() const noexcept { return _block_cache->hit_rate(); }
inline void HiCFile::reset_cache_stats() const noexcept { _block_cache->reset_stats(); }
inline void HiCFile::clear_cache() noexcept { _block_cache->clear(); }
inline void HiCFile::optimize_cache_size(std::size_t upper_bound) {
  return optimize_cache_size_for_random_access(upper_bound);
}

inline void HiCFile::optimize_cache_size_for_iteration(std::size_t upper_bound) {
  std::size_t cache_size = estimate_cache_size_cis() + estimate_cache_size_trans();
  cache_size = std::max(cache_size, std::size_t(10'000'000));
  _block_cache->set_capacity(std::min(upper_bound, cache_size));
}

inline void HiCFile::optimize_cache_size_for_random_access(std::size_t upper_bound) {
  std::size_t cache_size = estimate_cache_size_cis();
  cache_size = std::max(cache_size, std::size_t(10'000'000));
  _block_cache->set_capacity(std::min(upper_bound, cache_size));
}

inline std::size_t HiCFile::cache_capacity() const noexcept {
  return _block_cache->capacity_bytes();
}

inline std::size_t HiCFile::estimate_cache_size_cis() const {
  if (chromosomes().empty()) {
    return 0;
  }
  const auto& chrom1 = chromosomes().longest_chromosome();
  return fetch(chrom1.name(), chrom1.name()).estimate_optimal_cache_size();
}

inline std::size_t HiCFile::estimate_cache_size_trans() const {
  auto chrom1 = chromosomes().longest_chromosome();

  auto it = std::find_if(chromosomes().begin(), chromosomes().end(), [&](const Chromosome& chrom) {
    return !chrom.is_all() && chrom != chrom1;
  });
  if (it == chromosomes().end()) {
    return 0;
  }

  auto chrom2 = *it;

  if (chrom1.id() > chrom2.id()) {
    std::swap(chrom1, chrom2);
  }

  auto cache_size = fetch(chrom1.name(), chrom2.name()).estimate_optimal_cache_size();
  const auto num_trans_bins = bins().size() - bins().subset(chrom1).size();
  const auto num_chrom2_bins = bins().subset(chrom2).size();

  return ((cache_size + num_chrom2_bins - 1) / num_chrom2_bins) * num_trans_bins;
}

}  // namespace hictk::hic
