// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "hictk/hic/block_cache.hpp"
#include "hictk/hic/block_reader.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/hic/footer_cache.hpp"
#include "hictk/hic/header.hpp"
// #include "hictk/hic/hic_matrix_selector.hpp"
#include "hictk/hic/pixel_selector.hpp"

namespace hictk::hic {

class HiCFile {
  mutable std::shared_ptr<internal::HiCFileReader> _fs{};
  mutable internal::FooterCache _footers{};
  MatrixType _type{MatrixType::observed};
  MatrixUnit _unit{MatrixUnit::BP};
  mutable std::shared_ptr<internal::BlockCache> _block_cache{};
  std::shared_ptr<const BinTable> _bins{};

 public:
  using QUERY_TYPE = GenomicInterval::Type;
  explicit HiCFile(std::string url_, std::uint32_t resolution_,
                   MatrixType type_ = MatrixType::observed, MatrixUnit unit_ = MatrixUnit::BP,
                   std::uint64_t block_cache_capacity = 10'000'000);

  [[nodiscard]] HiCFile open_resolution(std::uint32_t resolution) const;
  [[nodiscard]] bool has_resolution(std::uint32_t resolution) const;

  [[nodiscard]] const std::string &url() const noexcept;
  [[nodiscard]] const std::string &name() const noexcept;
  [[nodiscard]] std::int32_t version() const noexcept;
  [[nodiscard]] const Reference &chromosomes() const noexcept;
  [[nodiscard]] const BinTable &bins() const noexcept;
  [[nodiscard]] const std::string &assembly() const noexcept;
  [[nodiscard]] const std::vector<std::uint32_t> &avail_resolutions() const noexcept;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  [[nodiscard]] PixelSelectorAll fetch(NormalizationMethod norm = NormalizationMethod::NONE) const;

  [[nodiscard]] PixelSelector fetch(std::string_view query,
                                    NormalizationMethod norm = NormalizationMethod::NONE,
                                    QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  [[nodiscard]] PixelSelector fetch(std::string_view chrom_name, std::uint32_t start,
                                    std::uint32_t end,
                                    NormalizationMethod norm = NormalizationMethod::NONE) const;
  [[nodiscard]] PixelSelector fetch(std::string_view range1, std::string_view range2,
                                    NormalizationMethod norm = NormalizationMethod::NONE,
                                    QUERY_TYPE query_type = QUERY_TYPE::UCSC) const;
  [[nodiscard]] PixelSelector fetch(std::string_view chrom1_name, std::uint32_t start1,
                                    std::uint32_t end1, std::string_view chrom2_name,
                                    std::uint32_t start2, std::uint32_t end2,
                                    NormalizationMethod norm = NormalizationMethod::NONE) const;
  [[nodiscard]] std::size_t num_cached_footers() const noexcept;
  void purge_footer_cache();

  [[nodiscard]] double block_cache_hit_rate() const noexcept;
  void reset_cache_stats() const noexcept;

 private:
  [[nodiscard]] std::shared_ptr<const internal::HiCFooter> get_footer(
      const Chromosome &chrom1, const Chromosome &chrom2, MatrixType matrix_type,
      NormalizationMethod norm, MatrixUnit unit, std::uint32_t resolution) const;

  [[nodiscard]] PixelSelector fetch(const Chromosome &chrom1, std::uint32_t start1,
                                    std::uint32_t end1, const Chromosome &chrom2,
                                    std::uint32_t start2, std::uint32_t end2,
                                    NormalizationMethod norm) const;
};

namespace utils {
[[nodiscard]] bool is_hic_file(const std::filesystem::path &path);
}  // namespace utils

}  // namespace hictk::hic

#include "../../hic_file_impl.hpp"
#include "../../hic_file_utils_impl.hpp"
