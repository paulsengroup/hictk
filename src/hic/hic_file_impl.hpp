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

#include "hictk/hic/common.hpp"
#include "hictk/hic/footer.hpp"

namespace hictk::hic {

inline HiCFile::HiCFile(std::string url_, std::uint32_t resolution_, MatrixType type_,
                        MatrixUnit unit_, std::uint64_t block_cache_capacity)
    : _fs(std::make_shared<internal::HiCFileReader>(std::move(url_))),
      _type(type_),
      _unit(unit_),
      _block_cache(std::make_shared<internal::BlockLRUCache>(block_cache_capacity)),
      _bins(std::make_shared<const BinTable>(chromosomes(), resolution_)) {
  assert(block_cache_capacity != 0);
  if (!has_resolution(resolution())) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("file {} does not have interactions for resolution {}"), url(), resolution()));
  }
}

inline HiCFile HiCFile::open_resolution(std::uint32_t resolution) const {
  return HiCFile(url(), resolution, _type, _unit);
}

inline bool HiCFile::has_resolution(std::uint32_t resolution) const {
  const auto match = std::find(avail_resolutions().begin(), avail_resolutions().end(), resolution);
  return match != avail_resolutions().end();
}

inline const std::string& HiCFile::url() const noexcept { return _fs->url(); }

inline const std::string& HiCFile::name() const noexcept { return url(); }

inline std::int32_t HiCFile::version() const noexcept { return _fs->version(); }

inline const Reference& HiCFile::chromosomes() const noexcept { return _fs->header().chromosomes; }

inline const std::string& HiCFile::assembly() const noexcept { return _fs->header().genomeID; }

inline const std::vector<std::uint32_t>& HiCFile::avail_resolutions() const noexcept {
  return _fs->header().resolutions;
}

inline std::uint32_t HiCFile::resolution() const noexcept { return _bins->bin_size(); }

inline std::shared_ptr<const internal::HiCFooter> HiCFile::get_footer(
    std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
    NormalizationMethod norm, MatrixUnit unit, std::uint32_t resolution) const {
  const internal::HiCFooterMetadata metadata{url(),
                                             matrix_type,
                                             norm,
                                             unit,
                                             resolution,
                                             _fs->header().chromosomes.at(chrom1_id),
                                             _fs->header().chromosomes.at(chrom2_id)};
  auto it = _footers.find(metadata);
  if (it != _footers.end()) {
    return *it;
  }

  auto [node, _] =
      _footers.emplace(_fs->read_footer(chrom1_id, chrom2_id, matrix_type, norm, unit, resolution));

  return *node;
}

inline PixelSelector HiCFile::fetch(std::string_view query, NormalizationMethod norm,
                                    QUERY_TYPE query_type,
                                    std::size_t read_at_once_threshold) const {
  const auto gi = query_type == QUERY_TYPE::BED
                      ? GenomicInterval::parse_bed(this->chromosomes(), query)
                      : GenomicInterval::parse_ucsc(this->chromosomes(), std::string{query});

  return this->fetch(gi.chrom(), gi.start(), gi.end(), gi.chrom(), gi.start(), gi.end(), norm,
                     read_at_once_threshold);
}

inline PixelSelector HiCFile::fetch(std::string_view chrom_name, std::uint32_t start,
                                    std::uint32_t end, NormalizationMethod norm,
                                    std::size_t read_at_once_threshold) const {
  return this->fetch(chrom_name, start, end, chrom_name, start, end, norm, read_at_once_threshold);
}

inline PixelSelector HiCFile::fetch(std::string_view range1, std::string_view range2,
                                    NormalizationMethod norm, QUERY_TYPE query_type,
                                    std::size_t read_at_once_threshold) const {
  const auto gi1 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(this->chromosomes(), range1)
                       : GenomicInterval::parse_ucsc(this->chromosomes(), std::string{range1});

  const auto gi2 = query_type == QUERY_TYPE::BED
                       ? GenomicInterval::parse_bed(this->chromosomes(), range2)
                       : GenomicInterval::parse_ucsc(this->chromosomes(), std::string{range2});

  return this->fetch(gi1.chrom(), gi1.start(), gi1.end(), gi2.chrom(), gi2.start(), gi2.end(), norm,
                     read_at_once_threshold);
}

inline PixelSelector HiCFile::fetch(std::string_view chrom1_name, std::uint32_t start1,
                                    std::uint32_t end1, std::string_view chrom2_name,
                                    std::uint32_t start2, std::uint32_t end2,
                                    NormalizationMethod norm,
                                    std::size_t read_at_once_threshold) const {
  return this->fetch(chromosomes().at(chrom1_name), start1, end1, chromosomes().at(chrom2_name),
                     start2, end2, norm, read_at_once_threshold);
}

inline PixelSelector HiCFile::fetch(const Chromosome& chrom1, std::uint32_t start1,
                                    std::uint32_t end1, const Chromosome& chrom2,
                                    std::uint32_t start2, std::uint32_t end2,
                                    NormalizationMethod norm,
                                    std::size_t read_at_once_threshold) const {
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

  auto footer = [&]() {
    try {
      return get_footer(chrom1.id(), chrom2.id(), _type, norm, _unit, resolution());
    } catch (const std::exception& e) {
      // Check whether query is valid but there are no interactions for the given chromosome
      // pair
      const auto missing_footer =
          std::string_view{e.what()}.find("unable to read file offset") == std::string_view::npos;
      if (missing_footer) {
        throw;
      }

      internal::HiCFooterMetadata metadata{url(),        _type,  norm,   _unit,
                                           resolution(), chrom1, chrom2, -1};

      if (metadata.fileOffset == -1) {
        return std::make_shared<const internal::HiCFooter>(internal::Index{metadata.chrom1,
                                                                           metadata.chrom2,
                                                                           metadata.unit,
                                                                           metadata.resolution,
                                                                           _fs->version(),
                                                                           1,
                                                                           1,
                                                                           0,
                                                                           {}},
                                                           std::move(metadata));
      }
      return std::make_shared<const internal::HiCFooter>(
          _fs->read_index(metadata.fileOffset, metadata.chrom1, metadata.chrom2, metadata.unit,
                          metadata.resolution),
          std::move(metadata));
    }
  }();

  return PixelSelector{_fs, footer, _block_cache, _bins, coord1, coord2, read_at_once_threshold};
}

inline std::size_t HiCFile::num_cached_footers() const noexcept { return _footers.size(); }

inline void HiCFile::purge_footer_cache() { _footers.clear(); }

inline double HiCFile::block_cache_hit_rate() const noexcept { return _block_cache->hit_rate(); }
inline std::size_t HiCFile::block_cache_size() const noexcept { return _block_cache->size(); }
inline void HiCFile::clear_block_cache() noexcept { _block_cache->reset(); }
}  // namespace hictk::hic
