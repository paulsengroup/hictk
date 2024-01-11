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
#include "hictk/hash.hpp"
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

struct HiCSectionOffsets {
  std::streamoff position;
  std::size_t size;
};

struct BlockIndexKey {
  Chromosome chrom1;
  Chromosome chrom2;
  std::uint32_t resolution;

  [[nodiscard]] bool operator<(const BlockIndexKey& other) const noexcept;
};

class MatrixBodyMetadataTank {
 public:
  struct Key {
    Chromosome chrom1{};
    Chromosome chrom2{};

    bool operator==(const Key& other) const noexcept;
  };

 private:
  phmap::flat_hash_map<Key, MatrixBodyMetadata> _tank{};
  phmap::flat_hash_map<Key, HiCSectionOffsets> _offsets{};

 public:
  MatrixBodyMetadataTank() = default;

  [[nodiscard]] bool contains(const Chromosome& chrom1, const Chromosome& chrom2) const noexcept;
  [[nodiscard]] auto at(const Chromosome& chrom1, const Chromosome& chrom2) const
      -> const MatrixBodyMetadata&;
  [[nodiscard]] HiCSectionOffsets offset(const Chromosome& chrom1, const Chromosome& chrom2) const;

  void insert(const Chromosome& chrom1, const Chromosome& chrom2, MatrixMetadata matrix_metadata,
              MatrixResolutionMetadata matrix_resolution_metadata);
  void update_offsets(const Chromosome& chrom1, const Chromosome& chrom2, std::streamoff position,
                      std::size_t size);

  void remove(const Chromosome& chrom1, const Chromosome& chrom2);

  auto operator()() const noexcept -> const phmap::flat_hash_map<Key, MatrixBodyMetadata>&;
  [[nodiscard]] std::string serialize(BinaryBuffer& buffer, bool clear = true) const;
};

class HiCFileWriter {
  std::shared_ptr<const HiCHeader> _header{};
  std::shared_ptr<filestream::FileStream> _fs{};
  std::unique_ptr<const hictk::internal::TmpDir> _tmpdir{};

  using BinTables = phmap::flat_hash_map<std::uint32_t, std::shared_ptr<const BinTable>>;
  using BlockIndex = phmap::btree_map<BlockIndexKey, phmap::btree_set<MatrixBlockMetadata>>;
  using BlockMappers = phmap::flat_hash_map<std::uint32_t, HiCInteractionToBlockMapper>;

  BinTables _bin_tables{};
  BlockIndex _block_index{};
  BlockMappers _block_mappers{};

  MatrixBodyMetadataTank _matrix_metadata{};
  using FooterTank = phmap::btree_map<std::pair<Chromosome, Chromosome>, FooterV5>;
  FooterTank _footers{};

  BinaryBuffer _bbuffer{};
  std::unique_ptr<libdeflate_compressor> _compressor{};
  std::string _compression_buffer{};

  HiCSectionOffsets _header_offsets{};
  HiCSectionOffsets _body_metadata_offsets{};
  HiCSectionOffsets _footer_offsets{};

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

  // Write header
  void write_header();
  void write_footer_offset();
  void write_norm_vector_index(std::streamoff position, std::size_t length);

  // Write pixels
  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void add_pixels(PixelIt first_pixel, PixelIt last_pixel);
  void write_pixels();

  // Write body
  auto write_body_metadata() -> HiCSectionOffsets;
  void add_body_metadata(std::uint32_t resolution, const Chromosome& chrom1,
                         const Chromosome& chrom2, const std::string& unit = "BP");

  // Write footer
  auto write_footers() -> HiCSectionOffsets;
  void add_footer(const Chromosome& chrom1, const Chromosome& chrom2);

  // Write expected/normalization values
  void write_expected_values(std::string_view unit);

  void finalize();

 private:
  [[nodiscard]] static auto init_bin_tables(const Reference& chromosomes,
                                            const std::vector<std::uint32_t>& resolutions)
      -> BinTables;
  [[nodiscard]] static auto init_interaction_block_mappers(const std::filesystem::path& root_folder,
                                                           const BinTables& bin_tables,
                                                           int compression_lvl) -> BlockMappers;
  template <typename PixelIt, typename = std::enable_if_t<is_iterable_v<PixelIt>>>
  void add_pixels(std::uint32_t resolution, PixelIt first_pixel, PixelIt last_pixel);

  auto write_pixels(const Chromosome& chrom1, const Chromosome& chrom2) -> HiCSectionOffsets;
  auto write_pixels(const Chromosome& chrom1, const Chromosome& chrom2, std::uint32_t resolution)
      -> HiCSectionOffsets;

  auto write_interaction_block(std::uint64_t block_id, const Chromosome& chrom1,
                               const Chromosome& chrom2, std::uint32_t resolution,
                               const MatrixInteractionBlock<float>& blk) -> HiCSectionOffsets;

  [[nodiscard]] std::size_t compute_block_column_count(const Chromosome& chrom1,
                                                       const Chromosome& chrom2,
                                                       std::uint32_t resolution);
  [[nodiscard]] std::size_t compute_num_bins(const Chromosome& chrom1, const Chromosome& chrom2,
                                             std::uint32_t resolution);
};
}  // namespace hictk::hic::internal

template <>
struct std::hash<hictk::hic::internal::MatrixBodyMetadataTank::Key> {
  inline std::size_t operator()(
      hictk::hic::internal::MatrixBodyMetadataTank::Key const& k) const noexcept {
    return hictk::internal::hash_combine(0, k.chrom1, k.chrom2);
  }
};

#include "./impl/file_writer_impl.hpp"  // NOLINT
