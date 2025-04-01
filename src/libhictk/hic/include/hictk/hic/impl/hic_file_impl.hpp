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
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/common.hpp"
#include "hictk/genomic_interval.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk::hic {

inline File::File(std::string url_, std::optional<std::uint32_t> resolution_, MatrixType type_,
                  MatrixUnit unit_, std::uint64_t block_cache_capacity)
    : _fs(std::make_shared<internal::HiCFileReader>(std::move(url_))),
      _type(type_),
      _unit(unit_),
      _block_cache(std::make_shared<internal::BlockCache>(block_cache_capacity)),
      _weight_cache(std::make_shared<internal::WeightCache>()),
      _bins(std::make_shared<const BinTable>(_fs->header().chromosomes,
                                             infer_or_validate_resolution(*_fs, resolution_))) {
  assert(has_resolution(resolution()));

  if (block_cache_capacity == 0) {
    optimize_cache_size();
  }
}

inline File& File::open(std::string url_, std::optional<std::uint32_t> resolution_,
                        MatrixType type_, MatrixUnit unit_, std::uint64_t block_cache_capacity) {
  if (_fs->path() == url_ && resolution() == resolution_.value_or(0) && _type == type_ &&
      _unit == unit_) {
    _block_cache->set_capacity(block_cache_capacity, false);
    return *this;
  }

  const auto prev_block_cache_capacity = _block_cache->capacity_bytes();
  *this = File(std::move(url_), resolution_, type_, unit_, block_cache_capacity);

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

inline std::uint32_t File::resolution() const noexcept { return bins().resolution(); }
inline std::uint64_t File::nbins() const { return bins().size(); }
inline std::uint64_t File::nchroms(bool include_ALL) const {
  if (include_ALL) {
    return chromosomes().size();
  }
  return chromosomes().remove_ALL().size();
}

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
  return _fs->list_avail_normalizations(_type, _unit, _bins->resolution());
}

inline std::shared_ptr<const internal::HiCFooter> File::get_footer(const Chromosome& chrom1,
                                                                   const Chromosome& chrom2,
                                                                   MatrixType matrix_type,
                                                                   const balancing::Method& norm,
                                                                   MatrixUnit unit) const {
  const internal::HiCFooterMetadata metadata{path(),       matrix_type, norm,  unit,
                                             resolution(), chrom1,      chrom2};
  auto it = _footers.find(metadata);
  if (it != _footers.end()) {
    return *it;
  }

  auto weights1 = _weight_cache->get_or_init(chrom1, norm);
  auto weights2 = _weight_cache->get_or_init(chrom2, norm);

  auto [node, _] = _footers.emplace(
      _fs->read_footer(chrom1, chrom2, *_bins, matrix_type, norm, unit, weights1, weights2));

  return *node;
}

constexpr auto File::matrix_type() const noexcept -> MatrixType { return _type; }

constexpr auto File::matrix_unit() const noexcept -> MatrixUnit { return _unit; }

inline PixelSelectorAll File::fetch(const balancing::Method& norm,
                                    std::optional<std::uint64_t> diagonal_band_width) const {
  std::vector<PixelSelector> selectors;
  bool file_is_empty = true;

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
        auto sel = fetch(chrom1.name(), chrom2.name(), norm, QUERY_TYPE::UCSC, diagonal_band_width);
        if (!sel.empty()) {
          selectors.emplace_back(std::move(sel));
          file_is_empty = false;
        }
      } catch (const std::exception& e) {
        const std::string_view msg{e.what()};
        const auto missing_norm = msg.find("unable to find") != std::string_view::npos &&
                                  msg.find("normalization vector") != std::string_view::npos;
        if (!missing_norm) {
          throw;
        }
        file_is_empty = false;
      }
    }
  }

  if (file_is_empty) {
    return PixelSelectorAll{{}, _weight_cache};
  }

  if (selectors.empty()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find {} normalization vectors at {} ({})"),
                    norm.to_string(), resolution(), _unit));
  }

  return PixelSelectorAll{std::move(selectors), _weight_cache};
}

inline PixelSelector File::fetch(std::string_view range, const balancing::Method& norm,
                                 QUERY_TYPE query_type,
                                 std::optional<std::uint64_t> diagonal_band_width) const {
  const auto gi = query_type == QUERY_TYPE::BED
                      ? GenomicInterval::parse_bed(chromosomes(), range)
                      : GenomicInterval::parse_ucsc(chromosomes(), std::string{range});

  return fetch(gi.chrom(), gi.start(), gi.end(), gi.chrom(), gi.start(), gi.end(), norm,
               diagonal_band_width);
}

inline PixelSelector File::fetch(std::string_view chrom_name, std::uint32_t start,
                                 std::uint32_t end, const balancing::Method& norm,
                                 std::optional<std::uint64_t> diagonal_band_width) const {
  return fetch(chrom_name, start, end, chrom_name, start, end, norm, diagonal_band_width);
}

inline PixelSelector File::fetch(std::string_view range1, std::string_view range2,
                                 const balancing::Method& norm, QUERY_TYPE query_type,
                                 std::optional<std::uint64_t> diagonal_band_width) const {
  const auto gi1 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(chromosomes(), range1)
                       : GenomicInterval::parse_ucsc(chromosomes(), std::string{range1});

  const auto gi2 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(chromosomes(), range2)
                       : GenomicInterval::parse_ucsc(chromosomes(), std::string{range2});

  return fetch(gi1.chrom(), gi1.start(), gi1.end(), gi2.chrom(), gi2.start(), gi2.end(), norm,
               diagonal_band_width);
}

inline PixelSelector File::fetch(std::string_view chrom1_name, std::uint32_t start1,
                                 std::uint32_t end1, std::string_view chrom2_name,
                                 std::uint32_t start2, std::uint32_t end2,
                                 const balancing::Method& norm,
                                 std::optional<std::uint64_t> diagonal_band_width) const {
  return fetch(chromosomes().at(chrom1_name), start1, end1, chromosomes().at(chrom2_name), start2,
               end2, norm, diagonal_band_width);
}

inline PixelSelector File::fetch(const Chromosome& chrom1, std::uint32_t start1, std::uint32_t end1,
                                 const Chromosome& chrom2, std::uint32_t start2, std::uint32_t end2,
                                 const balancing::Method& norm,
                                 std::optional<std::uint64_t> diagonal_band_width) const {
  if (chrom1 > chrom2) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("query {}:{}-{}; {}:{}-{}; overlaps with the lower-triangle of the matrix"),
        chrom1.name(), start1, end1, chrom2.name(), start2, end2));
  }

  const PixelCoordinates coord1 = {_bins->at(chrom1, start1), _bins->at(chrom1, end1 - 1)};
  const PixelCoordinates coord2 = {_bins->at(chrom2, start2), _bins->at(chrom2, end2 - 1)};

  return {_fs,
          get_footer(chrom1, chrom2, _type, norm, _unit),
          _block_cache,
          _bins,
          coord1,
          coord2,
          diagonal_band_width};
}

inline PixelSelector File::fetch(std::uint64_t first_bin, std::uint64_t last_bin,
                                 const balancing::Method& norm,
                                 std::optional<std::uint64_t> diagonal_band_width) const {
  return fetch(first_bin, last_bin, first_bin, last_bin, norm, diagonal_band_width);
}

inline PixelSelector File::fetch(std::uint64_t first_bin1, std::uint64_t last_bin1,
                                 std::uint64_t first_bin2, std::uint64_t last_bin2,
                                 const balancing::Method& norm,
                                 std::optional<std::uint64_t> diagonal_band_width) const {
  const PixelCoordinates coord1{bins().at(first_bin1), bins().at(last_bin1 - 1)};
  const PixelCoordinates coord2{bins().at(first_bin2), bins().at(last_bin2 - 1)};

  return fetch(coord1.bin1.chrom().name(), coord1.bin1.start(), coord1.bin2.end() - 1,
               coord2.bin1.chrom().name(), coord2.bin1.start(), coord2.bin2.end() - 1, norm,
               diagonal_band_width);
}

inline const balancing::Weights& File::normalization(const balancing::Method& norm,
                                                     const Chromosome& chrom) const {
  const auto w = normalization_ptr(norm, chrom);
  assert(w);
  return *w;
}
inline const balancing::Weights& File::normalization(std::string_view norm,
                                                     const Chromosome& chrom) const {
  const auto w = normalization_ptr(norm, chrom);
  assert(w);
  return *w;
}
inline const balancing::Weights& File::normalization(const balancing::Method& norm) const {
  const auto w = normalization_ptr(norm);
  assert(w);
  return *w;
}
inline const balancing::Weights& File::normalization(std::string_view norm) const {
  const auto w = normalization_ptr(norm);
  assert(w);
  return *w;
}

inline std::shared_ptr<const balancing::Weights> File::normalization_ptr(
    const balancing::Method& norm, const Chromosome& chrom) const {
  assert(_weight_cache);
  const auto expected_length = (chrom.size() + bins().resolution() - 1) / bins().resolution();

  try {
    // This takes care of populating the weight cache when appropriate
    const auto weight_size = fetch(chrom.name(), norm).weights1().size();
    if (weight_size != expected_length) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("{} normalization vector for {} appears to be corrupted: "
                                 "expected {} values, found {}"),
                      norm, chrom.name(), expected_length, weight_size));
    }
    return _weight_cache->at(chrom, norm);
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

  auto weights = _weight_cache->get_or_init(chrom, norm);
  assert(weights->empty());

  *weights = balancing::Weights{std::numeric_limits<double>::quiet_NaN(), expected_length,
                                balancing::Weights::Type::DIVISIVE};

  return weights;
}

inline std::shared_ptr<const balancing::Weights> File::normalization_ptr(
    std::string_view norm, const Chromosome& chrom) const {
  return normalization_ptr(balancing::Method{norm}, chrom);
}

inline std::shared_ptr<const balancing::Weights> File::normalization_ptr(
    const balancing::Method& norm) const {
  assert(_weight_cache);
  auto weights = _weight_cache->get_or_init(0, norm);
  if (!weights->empty()) {
    return weights;
  }

  if (norm == balancing::Method::NONE()) {
    *weights = balancing::Weights{1.0, bins().size(), balancing::Weights::Type::DIVISIVE};
    return weights;
  }

  std::vector<double> buff(bins().size(), std::numeric_limits<double>::quiet_NaN());
  for (const auto& chrom : chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }

    const auto& chrom_weights = normalization(norm, chrom);
    const auto offset = static_cast<std::ptrdiff_t>(bins().at(chrom).id());
    std::copy(chrom_weights.begin(balancing::Weights::Type::DIVISIVE),
              chrom_weights.end(balancing::Weights::Type::DIVISIVE), buff.begin() + offset);
  }

  *weights = balancing::Weights{std::move(buff), balancing::Weights::Type::DIVISIVE};
  return weights;
}

inline std::shared_ptr<const balancing::Weights> File::normalization_ptr(
    std::string_view norm) const {
  return normalization_ptr(balancing::Method{norm});
}

inline std::vector<double> File::expected_values(const Chromosome& chrom,
                                                 const balancing::Method& normalization_) const {
  const File f(path(), resolution(), MatrixType::expected, _unit, 1);
  const auto sel = f.fetch(chrom.name(), normalization_);
  const auto footer =
      f._footers.find(internal::HiCFooterMetadata{f.path(), MatrixType::expected, normalization_,
                                                  MatrixUnit::BP, resolution(), chrom, chrom, -1});
  if (footer == f._footers.end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to fetch expected values for \"{}\" ({})"), chrom.name(),
                    normalization_.to_string()));
  }
  return (*footer)->expectedValues();
}

inline std::size_t File::num_cached_footers() const noexcept { return _footers.size(); }

inline void File::purge_footer_cache() { _footers.clear(); }

inline double File::block_cache_hit_rate() const noexcept { return _block_cache->hit_rate(); }
inline void File::reset_cache_stats() const noexcept { _block_cache->reset_stats(); }
inline void File::clear_cache() noexcept { _block_cache->clear(); }
inline void File::optimize_cache_size(std::size_t upper_bound) {
  optimize_cache_size_for_random_access(upper_bound);
}

inline void File::optimize_cache_size_for_iteration(std::size_t upper_bound) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  if (version() < 9) {
    _block_cache->set_capacity(std::min(upper_bound, std::size_t{10'000'000}));
    return;
  }
  std::size_t cache_size = estimate_cache_size_cis() + estimate_cache_size_trans();
  cache_size = std::max(cache_size, std::size_t{10'000'000});
  _block_cache->set_capacity(std::min(upper_bound, cache_size));
  // NOLINTEND(*-avoid-magic-numbers)
}

inline void File::optimize_cache_size_for_random_access(std::size_t upper_bound) {
  // NOLINTBEGIN(*-avoid-magic-numbers)
  if (version() < 9) {
    _block_cache->set_capacity(std::min(upper_bound, std::size_t{10'000'000}));
    return;
  }
  std::size_t cache_size = estimate_cache_size_cis();
  cache_size = std::max(cache_size, std::size_t{10'000'000});
  _block_cache->set_capacity(std::min(upper_bound, cache_size));
  // NOLINTEND(*-avoid-magic-numbers)
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

inline std::uint32_t File::infer_or_validate_resolution(
    const internal::HiCFileReader& fs, std::optional<std::uint32_t> wanted_resolution) {
  const auto& resolutions = fs.header().resolutions;
  if (!wanted_resolution.has_value()) {
    if (resolutions.size() == 1) {
      return resolutions.front();
    }
    throw std::runtime_error("resolution is required when opening multi-resolution .hic files");
  }

  const auto match = std::find(resolutions.begin(), resolutions.end(), *wanted_resolution);

  if (match == resolutions.end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("file {} does not have interactions for resolution {}"), fs.path(),
                    *wanted_resolution));
  }

  return *wanted_resolution;
}

}  // namespace hictk::hic
