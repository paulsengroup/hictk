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

namespace hictk {

inline HiCFile::HiCFile(std::string url_, std::uint32_t resolution_, MatrixType type_,
                        MatrixUnit unit_, std::uint64_t block_cache_capacity)
    : _fs(std::make_shared<internal::HiCFileStream>(std::move(url_))),
      _type(type_),
      _unit(unit_),
      _block_cache(block_cache_capacity),
      _bins(chromosomes(), resolution_) {
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

constexpr std::uint32_t HiCFile::resolution() const noexcept { return _bins.bin_size(); }

inline std::shared_ptr<const internal::HiCFooter> HiCFile::get_footer(
    std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
    NormalizationMethod norm, MatrixUnit unit, std::uint32_t resolution) {
  const internal::HiCFooterMetadata metadata{url(),
                                             matrix_type,
                                             norm,
                                             unit,
                                             resolution,
                                             _fs->header().chromosomes.at(chrom1_id),
                                             _fs->header().chromosomes.at(chrom2_id)};
  auto it = _footers.find(metadata);
  if (it != _footers.end()) {
    return it->second;
  }
  auto footer = std::make_shared<const internal::HiCFooter>(
      _fs->readFooter(chrom1_id, chrom2_id, matrix_type, norm, unit, resolution));
  auto node = _footers.emplace(std::move(metadata), std::move(footer));

  assert(node.second);
  return node.first->second;
}

inline internal::MatrixSelector HiCFile::get_matrix_selector(const Chromosome& chrom,
                                                             NormalizationMethod norm) {
  return get_matrix_selector(chrom, chrom, norm);
}
inline internal::MatrixSelector HiCFile::get_matrix_selector(const std::string& chromName,
                                                             NormalizationMethod norm) {
  return get_matrix_selector(chromName, chromName, norm);
}
inline internal::MatrixSelector HiCFile::get_matrix_selector(std::uint32_t chrom_id,
                                                             NormalizationMethod norm) {
  return get_matrix_selector(chrom_id, chrom_id, norm);
}

inline internal::MatrixSelector HiCFile::get_matrix_selector(const Chromosome& chrom1,
                                                             const Chromosome& chrom2,
                                                             NormalizationMethod norm) {
  return get_matrix_selector(chrom1.id(), chrom2.id(), norm);
}

inline internal::MatrixSelector HiCFile::get_matrix_selector(const std::string& chrom1_name,
                                                             const std::string& chrom2_name,
                                                             NormalizationMethod norm) {
  const auto it1 = chromosomes().find(chrom1_name);
  if (it1 == chromosomes().end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find chromosome named {}"), chrom1_name));
  }
  if (chrom1_name == chrom2_name) {
    return get_matrix_selector(*it1, *it1, norm);
  }

  const auto it2 = chromosomes().find(chrom2_name);
  if (it2 == chromosomes().end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find chromosome named {}"), chrom2_name));
  }

  return get_matrix_selector(*it1, *it2, norm);
}

inline internal::MatrixSelector HiCFile::get_matrix_selector(std::uint32_t chrom1_id,
                                                             std::uint32_t chrom2_id,
                                                             NormalizationMethod norm) {
  if (chrom1_id >= std::int64_t(chromosomes().size())) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find chromosome corresponding to ID {}"), chrom1_id));
  }
  if (chrom2_id >= std::int64_t(chromosomes().size())) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find chromosome corresponding to ID {}"), chrom2_id));
  }

  if (chrom1_id > chrom2_id) {
    throw std::runtime_error(
        "Query overlaps the lower-triangle of the matrix. This is currently not supported.");
  }

  if (_type == MatrixType::expected && norm != NormalizationMethod::NONE) {
    throw std::logic_error(fmt::format(
        FMT_STRING("matrix type {} is incompatible with normalization method {}"), _type, norm));
  }

  try {
    return internal::MatrixSelector(
        _fs, get_footer(chrom1_id, chrom2_id, _type, norm, _unit, resolution()),
        1'000'000);  // TODO: REMOVE CACHE CAPACITY!
  } catch (const std::exception& e) {
    // Check whether query is valid but there are no interactions for the given chromosome pair
    const auto missing_footer =
        std::string_view{e.what()}.find("unable to read file offset") == std::string_view::npos;
    if (missing_footer) {
      throw;
    }

    internal::HiCFooterMetadata metadata{url(),
                                         _type,
                                         norm,
                                         _unit,
                                         resolution(),
                                         _fs->header().chromosomes.at(chrom1_id),
                                         _fs->header().chromosomes.at(chrom2_id),
                                         -1};

    return internal::MatrixSelector(
        _fs, std::make_shared<const internal::HiCFooter>(std::move(metadata)), 1);
  }
}

inline std::size_t HiCFile::num_cached_footers() const noexcept { return _footers.size(); }

inline void HiCFile::purge_footer_cache() { _footers.clear(); }

}  // namespace hictk
