// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <libdeflate.h>
#include <parallel_hashmap/btree.h>
#include <zstd.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <string>

#include "hictk/bin_table.hpp"
#include "hictk/hic/binary_buffer.hpp"
#include "hictk/hic/expected_values_aggregator.hpp"
#include "hictk/hic/file_writer_data_structures.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/hic/header.hpp"
#include "hictk/hic/interaction_block.hpp"
#include "hictk/hic/interaction_to_block_mapper.hpp"
#include "hictk/tmpdir.hpp"

template <>
struct std::default_delete<libdeflate_compressor> {
  void operator()(libdeflate_compressor* compressor) const {
    libdeflate_free_compressor(compressor);
  }
};

namespace hictk::hic::internal {

struct BlockIndexKey {
  Chromosome chrom1;
  Chromosome chrom2;
  std::uint32_t resolution;

  [[nodiscard]] bool operator<(const BlockIndexKey& other) const noexcept;
};

template <typename N>
class PixelTank {
  using ChromosomePair = std::pair<Chromosome, Chromosome>;
  using ChromPixelTank = std::vector<ThinPixel<N>>;

  std::shared_ptr<const BinTable> _bin_table{};
  phmap::flat_hash_map<ChromosomePair, ChromPixelTank> _pixel_tank{};
  phmap::flat_hash_map<ChromosomePair, float> _matrix_tot_counts{};
  ExpectedValuesAggregator _expected_values{};

 public:
  PixelTank() = default;
  explicit PixelTank(std::shared_ptr<const BinTable> bins);

  [[nodiscard]] bool contains(const Chromosome& chrom1, const Chromosome& chrom2) const noexcept;

  void add_pixel(const ThinPixel<N>& p, bool update_expected_values = true);
  void add_pixel(const Pixel<N>& p, bool update_expected_values = true);
  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void add_pixels(PixelIt first_pixel, PixelIt last_pixel, bool update_expected_values = true);

  void finalize();

  [[nodiscard]] auto pixels() const noexcept
      -> const phmap::flat_hash_map<ChromosomePair, ChromPixelTank>&;
  [[nodiscard]] auto matrix_counts() const noexcept
      -> const phmap::flat_hash_map<ChromosomePair, float>&;
  [[nodiscard]] auto expected_values() const noexcept;

  [[nodiscard]] auto pixels(const Chromosome& chrom1, const Chromosome& chrom2) const
      -> const ChromPixelTank&;
  [[nodiscard]] auto matrix_counts(const Chromosome& chrom1, const Chromosome& chrom2) const
      -> const float&;
};

class MetadataOffsetTank {
  struct Key {
    Chromosome chrom1{};
    Chromosome chrom2{};
    std::uint32_t resolution{};

    bool operator<(const Key& other) const noexcept;
  };

  struct Value {
    std::streamoff matrix_metadata_offset{};
    std::size_t matrix_size{};
  };

  phmap::btree_map<Key, Value> _tank{};

 public:
  MetadataOffsetTank() = default;

  [[nodiscard]] bool contains(const Chromosome& chrom1, const Chromosome& chrom2,
                              std::uint32_t resolution) const noexcept;
  [[nodiscard]] auto at(const Chromosome& chrom1, const Chromosome& chrom2,
                        std::uint32_t resolution) const -> const Value&;

  void insert(const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t resolution,
              std::size_t offset, std::size_t size);

  auto operator()() const noexcept -> const phmap::btree_map<Key, Value>&;
};

struct HiCSectionOffsets {
  std::streamoff position;
  std::size_t size;
};

class ChromChromHiCFileWriter {
  Chromosome _chrom1{};
  Chromosome _chrom2{};

  std::shared_ptr<const HiCHeader> _header{};
  std::shared_ptr<filestream::FileStream> _fs{};
  std::shared_ptr<const BinTable> _bin_table{};
  phmap::btree_map<BlockIndexKey, phmap::btree_set<MatrixBlockMetadata>> _block_index{};
  MatrixMetadata _matrix_metadata{};
  MatrixResolutionMetadata _matrix_resolution_metadata{};
  FooterV5 _footer{};

  BinaryBuffer _bbuffer{};
  std::unique_ptr<libdeflate_compressor> _compressor{};
  std::string _compression_buffer{};

  PixelTank<float> _pixel_tank{};

  std::unique_ptr<const hictk::internal::TmpDir> _tmpdir{};

  HiCSectionOffsets _header_offset{};
  HiCSectionOffsets _body_offset{};
  HiCSectionOffsets _body_metadata_offset{};

 public:
  ChromChromHiCFileWriter() = default;
  explicit ChromChromHiCFileWriter(
      const Chromosome& chrom1, const Chromosome& chrom2, HiCHeader header,
      const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
      std::int32_t compression_lvl = 9, std::size_t buffer_size = 32'000'000);

  [[nodiscard]] std::string_view url() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins() const;
  [[nodiscard]] std::uint32_t resolution() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void write_pixels(PixelIt first_pixel, PixelIt last_pixel, bool update_expected_values = true);

  // Write header
  void write_header();
  void write_footer_offset(std::streamoff master_index_offset);
  void write_norm_vector_index(std::streamoff position, std::size_t length);

  // Write body
  auto write_body_metadata(const std::string& unit = "BP") -> HiCSectionOffsets;

  auto write_interaction_block(std::uint64_t block_id, const std::vector<ThinPixel<float>>& pixels,
                               std::size_t bin_column_offset, std::size_t bin_row_offset)
      -> HiCSectionOffsets;

  // Write footer
  auto write_footer() -> HiCSectionOffsets;
  void write_footer_section_size(std::streamoff footer_offset, std::uint64_t bytes);

  // Write expected/normalization values
  void write_expected_values(std::string_view unit);

  // copy sections
  void copy_body(std::ofstream& dest, std::size_t chunk_size = 32'000'000);

  void finalize();

 private:
  [[nodiscard]] static std::shared_ptr<const HiCHeader> init_header(const Chromosome& chrom1,
                                                                    const Chromosome& chrom2,
                                                                    HiCHeader&& header);

  auto write_matrix_metadata() -> HiCSectionOffsets;
  auto write_resolution_metadata(float sum_counts, const std::string& unit) -> HiCSectionOffsets;
};

class HiCFileWriter {
  std::shared_ptr<const HiCHeader> _header{};
  std::shared_ptr<filestream::FileStream> _fs{};
  phmap::flat_hash_map<std::uint32_t, std::shared_ptr<const BinTable>> _bin_tables{};
  phmap::btree_map<BlockIndexKey, phmap::btree_set<MatrixBlockMetadata>> _block_index{};
  std::vector<MatrixMetadata> _matrix_metadata{};
  std::vector<MatrixResolutionMetadata> _matrix_resolution_metadata{};
  std::vector<FooterV5> _footers{};

  BinaryBuffer _bbuffer{};
  std::unique_ptr<libdeflate_compressor> _compressor{};
  std::string _compression_buffer{};

  MetadataOffsetTank _metadata_offset_tank{};
  phmap::flat_hash_map<std::uint32_t, PixelTank<float>> _pixel_tank{};

  std::unique_ptr<const hictk::internal::TmpDir> _tmpdir{};

  struct HiCSectionOffsets {
    std::streamoff position;
    std::size_t size;
  };

  HiCSectionOffsets _header_offsets{};
  HiCSectionOffsets _body_offsets{};
  [[maybe_unused]] HiCSectionOffsets _footer_offsets{};  // TODO

  static constexpr std::int32_t DEFAULT_INTRA_CUTOFF = 500;
  static constexpr std::int32_t DEFAULT_INTER_CUTOFF = 5'000;
  static constexpr std::size_t DEFAULT_BLOCK_CAPACITY = 1'000;

 public:
  HiCFileWriter() = default;
  explicit HiCFileWriter(
      HiCHeader header,
      const std::filesystem::path& tmpdir = std::filesystem::temp_directory_path(),
      std::int32_t compression_lvl = 9, std::size_t buffer_size = 32'000'000);

  [[nodiscard]] std::string_view url() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const BinTable& bins(std::uint32_t resolution) const;
  [[nodiscard]] const std::vector<std::uint32_t>& resolutions() const noexcept;

  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void append_pixels(std::uint32_t resolution, PixelIt first_pixel, PixelIt last_pixel,
                     bool update_expected_values = true);

  void write_pixels(bool write_chromosome_ALL = true);
  void write_pixels(const Chromosome& chrom1, const Chromosome& chrom2);
  void write_pixels(const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t resolution);
  void write_pixels_ALL(std::size_t num_bins = 1000);

  // Write header
  void write_header();
  void write_footer_offset(std::streamoff master_index_offset);
  void write_norm_vector_index(std::streamoff position, std::size_t length);

  // Write body
  auto write_body_metadata(std::uint32_t resolution, const Chromosome& chrom1,
                           const Chromosome& chrom2, const std::string& unit = "BP")
      -> HiCSectionOffsets;

  auto write_interaction_block(std::uint64_t block_id, const Chromosome& chrom1,
                               const Chromosome& chrom2, std::uint32_t resolution,
                               const std::vector<ThinPixel<float>>& pixels,
                               std::size_t bin_column_offset, std::size_t bin_row_offset)
      -> HiCSectionOffsets;

  // Write footer
  auto write_footers() -> HiCSectionOffsets;
  void add_footer(const Chromosome& chrom1, const Chromosome& chrom2, std::size_t file_offset,
                  std::size_t matrix_metadata_bytes);
  void write_footer_section_size(std::streamoff footer_offset, std::uint64_t bytes);

  // Write expected/normalization values
  void write_expected_values(std::string_view unit);

  // copy sections
  void copy_body(HiCFileWriter& dest, std::size_t chunk_size = 32'000'000);

  void finalize();

 private:
  [[nodiscard]] std::size_t compute_block_column_count(
      std::size_t num_bins, std::uint32_t bin_size, std::uint32_t cutoff,
      std::size_t block_capacity = DEFAULT_BLOCK_CAPACITY);
  [[nodiscard]] std::size_t compute_num_bins(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                             std::size_t bin_size);

  auto write_matrix_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id) -> HiCSectionOffsets;
  auto write_resolutions_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                  float sum_counts, const std::string& unit) -> HiCSectionOffsets;

  void coarsen_pixels();
};
}  // namespace hictk::hic::internal

#include "./impl/file_writer_impl.hpp"  // NOLINT
