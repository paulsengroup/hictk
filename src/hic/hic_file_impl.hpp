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

inline HiCFile::HiCFile(std::string url_)
    : _fs(std::make_shared<internal::HiCFileStream>(std::move(url_))) {}

inline const std::string& HiCFile::url() const noexcept { return _fs->url(); }

inline const std::string& HiCFile::name() const noexcept { return url(); }

inline std::int32_t HiCFile::version() const noexcept { return _fs->version(); }

inline const Reference& HiCFile::chromosomes() const noexcept { return _fs->header().chromosomes; }

inline const std::string& HiCFile::assembly() const noexcept { return _fs->header().genomeID; }

inline const std::vector<std::uint32_t>& HiCFile::resolutions() const noexcept {
  return _fs->header().resolutions;
}

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

inline internal::MatrixSelector HiCFile::get_matrix_selector(
    const Chromosome& chrom, MatrixType matrix_type, NormalizationMethod norm, MatrixUnit unit,
    std::uint32_t resolution, std::size_t block_cache_capacity) {
  return get_matrix_selector(chrom, chrom, matrix_type, norm, unit, resolution,
                             block_cache_capacity);
}
inline internal::MatrixSelector HiCFile::get_matrix_selector(
    const std::string& chromName, MatrixType matrix_type, NormalizationMethod norm, MatrixUnit unit,
    std::uint32_t resolution, std::size_t block_cache_capacity) {
  return get_matrix_selector(chromName, chromName, matrix_type, norm, unit, resolution,
                             block_cache_capacity);
}
inline internal::MatrixSelector HiCFile::get_matrix_selector(
    std::uint32_t chrom_id, MatrixType matrix_type, NormalizationMethod norm, MatrixUnit unit,
    std::uint32_t resolution, std::size_t block_cache_capacity) {
  return get_matrix_selector(chrom_id, chrom_id, matrix_type, norm, unit, resolution,
                             block_cache_capacity);
}

inline internal::MatrixSelector HiCFile::get_matrix_selector(
    const Chromosome& chrom1, const Chromosome& chrom2, MatrixType matrix_type,
    NormalizationMethod norm, MatrixUnit unit, std::uint32_t resolution,
    std::size_t block_cache_capacity) {
  return get_matrix_selector(chrom1.id(), chrom2.id(), matrix_type, norm, unit, resolution,
                             block_cache_capacity);
}

inline internal::MatrixSelector HiCFile::get_matrix_selector(
    const std::string& chrom1_name, const std::string& chrom2_name, MatrixType matrix_type,
    NormalizationMethod norm, MatrixUnit unit, std::uint32_t resolution,
    std::size_t block_cache_capacity) {
  const auto it1 = chromosomes().find(chrom1_name);
  if (it1 == chromosomes().end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find chromosome named {}"), chrom1_name));
  }
  if (chrom1_name == chrom2_name) {
    return get_matrix_selector(*it1, *it1, matrix_type, norm, unit, resolution,
                               block_cache_capacity);
  }

  const auto it2 = chromosomes().find(chrom2_name);
  if (it2 == chromosomes().end()) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("unable to find chromosome named {}"), chrom2_name));
  }

  return get_matrix_selector(*it1, *it2, matrix_type, norm, unit, resolution, block_cache_capacity);
}

inline internal::MatrixSelector HiCFile::get_matrix_selector(
    std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
    NormalizationMethod norm, MatrixUnit unit, std::uint32_t resolution,
    std::size_t block_cache_capacity) {
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

  if (matrix_type == MatrixType::expected && norm != NormalizationMethod::NONE) {
    throw std::logic_error(
        fmt::format(FMT_STRING("matrix type {} is incompatible with normalization method {}"),
                    matrix_type, norm));
  }

  const auto it = std::find(resolutions().begin(), resolutions().end(), resolution);
  if (it == resolutions().end()) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "matrix does not have interactions for resolution {}. Available resolutions: {}"),
        resolution, fmt::join(_fs->header().resolutions, ", ")));
  }

  try {
    return internal::MatrixSelector(
        _fs, get_footer(chrom1_id, chrom2_id, matrix_type, norm, unit, resolution),
        block_cache_capacity);
  } catch (const std::exception& e) {
    // Check whether query is valid but there are no interactions for the given chromosome pair
    const auto missing_footer =
        std::string_view{e.what()}.find("unable to read file offset") == std::string_view::npos;
    if (missing_footer) {
      throw;
    }

    internal::HiCFooterMetadata metadata{url(),
                                         matrix_type,
                                         norm,
                                         unit,
                                         resolution,
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
