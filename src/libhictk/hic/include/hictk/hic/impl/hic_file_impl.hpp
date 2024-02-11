// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk::hic {

inline File::File(std::string url_, std::uint32_t resolution_, MatrixType type_, MatrixUnit unit_,
                  std::uint64_t block_cache_capacity)
    : _fs(std::make_shared<internal::HiCFileReader>(std::move(url_))),
      _type(type_),
      _unit(unit_),
      _block_cache(std::make_shared<internal::BlockCache>(block_cache_capacity)),
      _weight_cache(std::make_shared<internal::WeightCache>()),
      _bins(std::make_shared<const BinTable>(_fs->header().chromosomes, resolution_)) {
  if (!has_resolution(resolution())) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("file {} does not have interactions for resolution {}"), path(), resolution()));
  }

  if (block_cache_capacity == 0) {
    optimize_cache_size();
  }
}

inline File& File::open(std::string url_, std::uint32_t resolution_, MatrixType type_,
                        MatrixUnit unit_, std::uint64_t block_cache_capacity) {
  if (_fs->path() == url_ && resolution() == resolution_ && _type == type_ && _unit == unit_) {
    _block_cache->set_capacity(block_cache_capacity, false);
    return *this;
  }

  const auto prev_block_cache_capacity = _block_cache->capacity_bytes();
  *this = File(url_, resolution_, type_, unit_, block_cache_capacity);

  if (_block_cache->capacity_bytes() < prev_block_cache_capacity) {
    _block_cache->set_capacity(prev_block_cache_capacity);
  }
  return *this;
}

inline File& File::open(std::uint32_t resolution_, MatrixType type_, MatrixUnit unit_,
                        std::uint64_t block_cache_capacity) {
  return open(path(), resolution_, type_, unit_, block_cache_capacity);
}

inline bool File::has_resolution(std::uint32_t resolution) const {
  const auto match = std::find(avail_resolutions().begin(), avail_resolutions().end(), resolution);
  return match != avail_resolutions().end();
}

inline const std::string& File::path() const noexcept { return _fs->path(); }

inline const std::string& File::name() const noexcept { return path(); }

inline std::int32_t File::version() const noexcept { return _fs->version(); }

inline const BinTable& File::bins() const noexcept {
  assert(_bins);
  return *_bins;
}

inline std::shared_ptr<const BinTable> File::bins_ptr() const noexcept { return _bins; }

inline std::uint32_t File::bin_size() const noexcept { return bins().bin_size(); }
inline std::uint64_t File::nbins() const { return bins().size(); }
inline std::uint64_t File::nchroms() const { return chromosomes().size(); }

inline const Reference& File::chromosomes() const noexcept { return bins().chromosomes(); }

inline const std::string& File::assembly() const noexcept { return _fs->header().genomeID; }

inline const phmap::flat_hash_map<std::string, std::string>& File::attributes() const noexcept {
  return _fs->header().attributes;
}

inline const std::vector<std::uint32_t>& File::avail_resolutions() const noexcept {
  return _fs->header().resolutions;
}

inline bool File::has_normalization(std::string_view normalization) const {
  const auto normalizations = avail_normalizations();
  const auto it = std::find_if(normalizations.begin(), normalizations.end(),
                               [&](const auto& norm) { return norm.to_string() == normalization; });

  return it != normalizations.end();
}

inline std::vector<balancing::Method> File::avail_normalizations() const {
  return _fs->list_avail_normalizations(_type, _unit, _bins->bin_size());
}

inline std::uint32_t File::resolution() const noexcept { return _bins->bin_size(); }

inline std::shared_ptr<const internal::HiCFooter> File::get_footer(
    const Chromosome& chrom1, const Chromosome& chrom2, MatrixType matrix_type,
    balancing::Method norm, MatrixUnit unit, std::uint32_t resolution) const {
  const internal::HiCFooterMetadata metadata{path(),     matrix_type, norm,  unit,
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

inline PixelSelectorAll File::fetch(balancing::Method norm) const {
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
      try {
        auto sel = fetch(chrom1.name(), chrom2.name(), norm);
        if (!sel.empty()) {
          selectors.emplace_back(std::move(sel));
        }
      } catch (const std::exception& e) {
        const std::string_view msg{e.what()};
        const auto missing_norm = msg.find("unable to find") != std::string_view::npos &&
                                  msg.find("normalization vector") != std::string_view::npos;
        if (!missing_norm) {
          throw;
        }
      }
    }
  }

  if (selectors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find {} normalization vectors at {} ({})"),
                    norm.to_string(), resolution(), _unit));
  }

  return PixelSelectorAll{std::move(selectors)};
}

inline PixelSelector File::fetch(std::string_view range, balancing::Method norm,
                                 QUERY_TYPE query_type) const {
  const auto gi = query_type == QUERY_TYPE::BED
                      ? GenomicInterval::parse_bed(chromosomes(), range)
                      : GenomicInterval::parse_ucsc(chromosomes(), std::string{range});

  return fetch(gi.chrom(), gi.start(), gi.end(), gi.chrom(), gi.start(), gi.end(), norm);
}

inline PixelSelector File::fetch(std::string_view chrom_name, std::uint32_t start,
                                 std::uint32_t end, balancing::Method norm) const {
  return fetch(chrom_name, start, end, chrom_name, start, end, norm);
}

inline PixelSelector File::fetch(std::string_view range1, std::string_view range2,
                                 balancing::Method norm, QUERY_TYPE query_type) const {
  const auto gi1 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(chromosomes(), range1)
                       : GenomicInterval::parse_ucsc(chromosomes(), std::string{range1});

  const auto gi2 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(chromosomes(), range2)
                       : GenomicInterval::parse_ucsc(chromosomes(), std::string{range2});

  return fetch(gi1.chrom(), gi1.start(), gi1.end(), gi2.chrom(), gi2.start(), gi2.end(), norm);
}

inline PixelSelector File::fetch(std::string_view chrom1_name, std::uint32_t start1,
                                 std::uint32_t end1, std::string_view chrom2_name,
                                 std::uint32_t start2, std::uint32_t end2,
                                 balancing::Method norm) const {
  return fetch(chromosomes().at(chrom1_name), start1, end1, chromosomes().at(chrom2_name), start2,
               end2, norm);
}

inline PixelSelector File::fetch(const Chromosome& chrom1, std::uint32_t start1, std::uint32_t end1,
                                 const Chromosome& chrom2, std::uint32_t start2, std::uint32_t end2,
                                 balancing::Method norm) const {
  if (chrom1 > chrom2) {
    throw std::runtime_error(
        "Query overlaps the lower-triangle of the matrix. This is currently not supported.");
  }

  const PixelCoordinates coord1 = {_bins->at(chrom1, start1), _bins->at(chrom1, end1 - 1)};
  const PixelCoordinates coord2 = {_bins->at(chrom2, start2), _bins->at(chrom2, end2 - 1)};

  return {_fs,          get_footer(chrom1, chrom2, _type, norm, _unit, resolution()),
          _block_cache, _bins,
          coord1,       coord2};
}

inline PixelSelector File::fetch(std::uint64_t first_bin, std::uint64_t last_bin,
                                 balancing::Method norm) const {
  return fetch(first_bin, last_bin, first_bin, last_bin, std::move(norm));
}

inline PixelSelector File::fetch(std::uint64_t first_bin1, std::uint64_t last_bin1,
                                 std::uint64_t first_bin2, std::uint64_t last_bin2,
                                 balancing::Method norm) const {
  PixelCoordinates coord1{bins().at(first_bin1), bins().at(last_bin1 - 1)};
  PixelCoordinates coord2{bins().at(first_bin2), bins().at(last_bin2 - 1)};

  return fetch(coord1.bin1.chrom().name(), coord1.bin1.start(), coord1.bin2.end() - 1,
               coord2.bin1.chrom().name(), coord2.bin1.start(), coord2.bin2.end() - 1,
               std::move(norm));
}

inline balancing::Weights File::normalization(balancing::Method norm,
                                              const Chromosome& chrom) const {
  std::vector<double> weights_{};
  const auto expected_length = (chrom.size() + bins().bin_size() - 1) / bins().bin_size();
  try {
    auto weights = fetch(chrom.name(), norm).weights1();
    if (!!weights && weights().size() != expected_length) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("{} normalization vector for {} appears to be corrupted: "
                                 "expected {} values, found {}"),
                      norm, chrom.name(), expected_length, weights().size()));
    }
    weights_ = weights();
  } catch (const std::exception& e) {
    const std::string_view msg{e.what()};

    const auto missing_interactions =
        msg.find("unable to read file offset") != std::string_view::npos;

    const auto missing_norm_vect =
        msg.find(fmt::format(FMT_STRING("unable to find {} normalization vector"), norm)) !=
        std::string_view::npos;

    if (!missing_interactions && !missing_norm_vect) {
      throw;
    }
  }

  if (weights_.empty()) {
    weights_.resize(expected_length, std::numeric_limits<double>::quiet_NaN());
  }

  return {weights_, balancing::Weights::Type::DIVISIVE};
}

inline balancing::Weights File::normalization(std::string_view norm,
                                              const Chromosome& chrom) const {
  return normalization(balancing::Method{norm}, chrom);
}

inline balancing::Weights File::normalization(balancing::Method norm) const {
  std::vector<double> weights{};
  weights.reserve(bins().size());
  for (const auto& chrom : chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }

    const auto chrom_weights = normalization(norm, chrom);
    weights.insert(weights.end(), chrom_weights().begin(), chrom_weights().end());
  }

  assert(weights.size() == bins().size());
  return {weights, balancing::Weights::Type::DIVISIVE};
}

inline balancing::Weights File::normalization(std::string_view norm) const {
  return normalization(balancing::Method{norm});
}

inline std::size_t File::num_cached_footers() const noexcept { return _footers.size(); }

inline void File::purge_footer_cache() { _footers.clear(); }

inline double File::block_cache_hit_rate() const noexcept { return _block_cache->hit_rate(); }
inline void File::reset_cache_stats() const noexcept { _block_cache->reset_stats(); }
inline void File::clear_cache() noexcept { _block_cache->clear(); }
inline void File::optimize_cache_size(std::size_t upper_bound) {
  return optimize_cache_size_for_random_access(upper_bound);
}

inline void File::optimize_cache_size_for_iteration(std::size_t upper_bound) {
  if (version() < 9) {
    _block_cache->set_capacity(std::min(upper_bound, std::size_t(10'000'000)));
    return;
  }
  std::size_t cache_size = estimate_cache_size_cis() + estimate_cache_size_trans();
  cache_size = std::max(cache_size, std::size_t(10'000'000));
  _block_cache->set_capacity(std::min(upper_bound, cache_size));
}

inline void File::optimize_cache_size_for_random_access(std::size_t upper_bound) {
  if (version() < 9) {
    _block_cache->set_capacity(std::min(upper_bound, std::size_t(10'000'000)));
    return;
  }
  std::size_t cache_size = estimate_cache_size_cis();
  cache_size = std::max(cache_size, std::size_t(10'000'000));
  _block_cache->set_capacity(std::min(upper_bound, cache_size));
}

inline std::size_t File::cache_capacity() const noexcept { return _block_cache->capacity_bytes(); }

inline std::size_t File::estimate_cache_size_cis() const {
  if (chromosomes().empty()) {
    return 0;
  }
  const auto& chrom1 = chromosomes().longest_chromosome();
  return fetch(chrom1.name(), chrom1.name()).estimate_optimal_cache_size();
}

inline std::size_t File::estimate_cache_size_trans() const {
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
