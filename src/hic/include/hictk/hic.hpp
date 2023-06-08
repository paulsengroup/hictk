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

 public:
  static constexpr std::size_t DEFAULT_BLOCK_CACHE_CAPACITY = 500UL << 20U;  // ~500MB
  explicit HiCFile(std::string url_);

  [[nodiscard]] const std::string &url() const noexcept;
  [[nodiscard]] const std::string &name() const noexcept;
  [[nodiscard]] std::int32_t version() const noexcept;
  [[nodiscard]] const ChromosomeMap &chromosomes() const noexcept;
  [[nodiscard]] const std::string &assembly() const noexcept;
  [[nodiscard]] const std::vector<std::int32_t> &resolutions() const noexcept;

  [[nodiscard]] internal::MatrixSelector get_matrix_selector(
      const chromosome &chrom, MatrixType matrix_type, NormalizationMethod norm, MatrixUnit unit,
      std::int32_t resolution, std::size_t block_cache_capacity = DEFAULT_BLOCK_CACHE_CAPACITY);
  [[nodiscard]] internal::MatrixSelector get_matrix_selector(
      const std::string &chromName, MatrixType matrix_type, NormalizationMethod norm,
      MatrixUnit unit, std::int32_t resolution,
      std::size_t block_cache_capacity = DEFAULT_BLOCK_CACHE_CAPACITY);
  [[nodiscard]] internal::MatrixSelector get_matrix_selector(
      std::int32_t chrom_id, MatrixType matrix_type, NormalizationMethod norm, MatrixUnit unit,
      std::int32_t resolution, std::size_t block_cache_capacity = DEFAULT_BLOCK_CACHE_CAPACITY);

  [[nodiscard]] internal::MatrixSelector get_matrix_selector(
      const chromosome &chrom1, const chromosome &chrom2, MatrixType matrix_type,
      NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution,
      std::size_t block_cache_capacity = DEFAULT_BLOCK_CACHE_CAPACITY);
  [[nodiscard]] internal::MatrixSelector get_matrix_selector(
      const std::string &chrom1_name, const std::string &chrom2_name, MatrixType matrix_type,
      NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution,
      std::size_t block_cache_capacity = DEFAULT_BLOCK_CACHE_CAPACITY);
  [[nodiscard]] internal::MatrixSelector get_matrix_selector(
      std::int32_t chrom1_id, std::int32_t chrom2_id, MatrixType matrix_type,
      NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution,
      std::size_t block_cache_capacity = DEFAULT_BLOCK_CACHE_CAPACITY);

  [[nodiscard]] std::size_t num_cached_footers() const noexcept;
  void purge_footer_cache();

 private:
  [[nodiscard]] std::shared_ptr<const internal::HiCFooter> get_footer(
      std::int32_t chrom1_id, std::int32_t chrom2_id, MatrixType matrix_type,
      NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution);
};

namespace utils {
[[nodiscard]] bool is_hic_file(const std::filesystem::path &path);
}  // namespace utils

}  // namespace hictk

#include "../../hic_file_impl.hpp"
#include "../../hic_file_utils_impl.hpp"
