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

#include "hictk/hic/common.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/hic_file_stream.hpp"
#include "hictk/hic/hic_footer.hpp"
#include "hictk/hic/hic_header.hpp"
#include "hictk/hic/hic_matrix_selector.hpp"

namespace hictk {
class HiCFile {
  // clang-format off
    using FooterCacheT =
        std::unordered_map<internal::HiCFooterMetadata,
                           std::shared_ptr<const internal::HiCFooter>>;
  // clang-format on
  std::shared_ptr<internal::HiCFileStream> _fs{};
  FooterCacheT _footers{};
  MatrixType _type{MatrixType::observed};
  MatrixUnit _unit{MatrixUnit::BP};
  internal::BlockLRUCache _block_cache{};
  BinTable _bins{};

 public:
  explicit HiCFile(std::string url_, std::uint32_t resolution_,
                   MatrixType type_ = MatrixType::observed, MatrixUnit unit_ = MatrixUnit::BP,
                   // TODO consider expressing cache size in terms of number of pixels
                   std::uint64_t block_cache_capacity = 500ULL << 20U);

  [[nodiscard]] HiCFile open_resolution(std::uint32_t resolution) const;
  [[nodiscard]] bool has_resolution(std::uint32_t resolution) const;

  [[nodiscard]] const std::string &url() const noexcept;
  [[nodiscard]] const std::string &name() const noexcept;
  [[nodiscard]] std::int32_t version() const noexcept;
  [[nodiscard]] const Reference &chromosomes() const noexcept;
  [[nodiscard]] const std::string &assembly() const noexcept;
  [[nodiscard]] const std::vector<std::uint32_t> &avail_resolutions() const noexcept;
  [[nodiscard]] constexpr std::uint32_t resolution() const noexcept;

  [[nodiscard]] internal::MatrixSelector get_matrix_selector(const Chromosome &chrom,
                                                             NormalizationMethod norm);
  [[nodiscard]] internal::MatrixSelector get_matrix_selector(const std::string &chromName,
                                                             NormalizationMethod norm);
  [[nodiscard]] internal::MatrixSelector get_matrix_selector(std::uint32_t chrom_id,
                                                             NormalizationMethod norm);

  [[nodiscard]] internal::MatrixSelector get_matrix_selector(const Chromosome &chrom1,
                                                             const Chromosome &chrom2,
                                                             NormalizationMethod norm);
  [[nodiscard]] internal::MatrixSelector get_matrix_selector(const std::string &chrom1_name,
                                                             const std::string &chrom2_name,
                                                             NormalizationMethod norm);
  [[nodiscard]] internal::MatrixSelector get_matrix_selector(std::uint32_t chrom1_id,
                                                             std::uint32_t chrom2_id,
                                                             NormalizationMethod norm);

  [[nodiscard]] std::size_t num_cached_footers() const noexcept;
  void purge_footer_cache();

 private:
  [[nodiscard]] std::shared_ptr<const internal::HiCFooter> get_footer(
      std::uint32_t chrom1_id, std::uint32_t chrom2_id, MatrixType matrix_type,
      NormalizationMethod norm, MatrixUnit unit, std::uint32_t resolution);
};

namespace utils {
[[nodiscard]] bool is_hic_file(const std::filesystem::path &path);
}  // namespace utils

}  // namespace hictk

#include "../../hic_file_impl.hpp"
#include "../../hic_file_utils_impl.hpp"
