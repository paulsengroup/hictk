// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/compile.h>
#include <fmt/format.h>
#include <libdeflate.h>
#include <parallel_hashmap/phmap.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <ios>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/reference.hpp"

namespace hictk::hic::internal {

inline std::string MatrixMetadata::serialize(BinaryBuffer &buffer) const {
  buffer.clear();
  buffer.write(chr1Idx);
  buffer.write(chr2Idx);
  buffer.write(nResolutions);

  return buffer.get();
}

inline std::string MatrixBlockMetadata::serialize(BinaryBuffer &buffer) const {
  buffer.clear();
  buffer.write(blockNumber);
  buffer.write(blockPosition);
  buffer.write(blockSizeBytes);

  return buffer.get();
}

inline bool MatrixBlockMetadata::operator<(const MatrixBlockMetadata &other) const noexcept {
  return blockNumber < other.blockNumber;
}

inline std::string MatrixResolutionMetadata::serialize(BinaryBuffer &buffer) const {
  buffer.clear();

  buffer.write(unit);
  buffer.write(resIdx);
  buffer.write(sumCounts);
  buffer.write(occupiedCellCount);
  buffer.write(percent5);
  buffer.write(percent95);
  buffer.write(binSize);
  buffer.write(blockSize);
  buffer.write(blockColumnCount);
  buffer.write(blockCount);

  return buffer.get();
}

inline MatrixInteractionBlock::MatrixInteractionBlock(const BinTable &bins,
                                                      const std::vector<ThinPixel<float>> &pixels,
                                                      std::size_t bin_row_offset)
    : nRecords(static_cast<std::int32_t>(pixels.size())),
      binRowOffset(static_cast<std::int32_t>(bin_row_offset)),
      _interactions(group_interactions_by_column(bins, pixels)) {}

inline std::string MatrixInteractionBlock::serialize(BinaryBuffer &buffer,
                                                     libdeflate_compressor &compressor,
                                                     std::string &compression_buffer) const {
  // TODO support dense layout
  // TODO support representation using shorts

  buffer.clear();

  buffer.write(nRecords);
  buffer.write(binColumnOffset);
  buffer.write(binRowOffset);
  buffer.write(useFloatContact);
  buffer.write(useIntXPos);
  buffer.write(useIntYPos);
  buffer.write(matrixRepresentation);

  const auto rowCount = static_cast<std::int32_t>(_interactions.size());  // TODO support short
  buffer.write(rowCount);

  for (const auto &[row, pixels] : _interactions) {
    assert(static_cast<std::int32_t>(row) >= binRowOffset);
    const auto rowNumber = static_cast<std::int32_t>(row) - binRowOffset;  // TODO support short
    const auto recordCount = static_cast<std::int32_t>(pixels.size());     // TODO support short
    buffer.write(rowNumber);
    buffer.write(recordCount);

    assert(std::is_sorted(pixels.begin(), pixels.end()));
    for (const auto &p : pixels) {
      const auto bin_id = static_cast<std::int32_t>(p.coords.bin1.rel_id());
      assert(bin_id >= binColumnOffset);
      const auto binColumn = bin_id - binColumnOffset;
      const auto value = p.count;
      buffer.write(binColumn);
      buffer.write(value);
    }
  }

  assert(compression_buffer.capacity() != 0);
  compression_buffer.resize(compression_buffer.capacity());
  while (true) {
    const auto compressed_size =
        libdeflate_zlib_compress(&compressor, buffer.get().data(), buffer.get().size(),
                                 compression_buffer.data(), compression_buffer.size());
    if (compressed_size != 0) {
      compression_buffer.resize(compressed_size);
      break;
    }

    compression_buffer.resize(compression_buffer.size() * 2);
  }

  return compression_buffer;
}

inline auto MatrixInteractionBlock::group_interactions_by_column(
    const BinTable &bins, const std::vector<ThinPixel<float>> &pixels)
    -> phmap::btree_map<RowID, Row> {
  phmap::btree_map<RowID, Row> buffer;

  for (const auto &p : pixels) {
    Pixel<float> pp(bins, p);
    const auto col = static_cast<std::int32_t>(pp.coords.bin2.rel_id());
    auto it = buffer.find(col);
    if (it != buffer.end()) {
      it->second.push_back(pp);
    } else {
      buffer.emplace(col, std::vector<Pixel<float>>{pp});
    }
  }

  for (auto &[_, v] : buffer) {
    std::sort(v.begin(), v.end());
  }
  return buffer;
}

inline std::string MasterIndex::serialize(BinaryBuffer &buffer) const {
  buffer.clear();
  buffer.write(key);
  buffer.write(position);
  buffer.write(size);

  return buffer.get();
}

inline std::string ExpectedValues::serialize(BinaryBuffer &buffer) const {
  buffer.clear();
  buffer.write(nExpectedValueVectors);

  return buffer.get();
}

inline std::string NormalizedExpectedValues::serialize(BinaryBuffer &buffer) const {
  buffer.clear();
  buffer.write(nNormExpectedValueVectors);

  return buffer.get();
}

inline std::string NormalizationVectorIndex::serialize(BinaryBuffer &buffer) const {
  buffer.clear();
  buffer.write(nNormVectors);

  return buffer.get();
}

inline std::string NormalizationVectorArray::serialize(BinaryBuffer &buffer) const {
  buffer.clear();
  buffer.write(nValues);

  return buffer.get();
}

inline std::string FooterV5::serialize(BinaryBuffer &buffer) const {
  std::string data = masterIndex.serialize(buffer);
  data += expectedValues.serialize(buffer);
  data += normExpectedValues.serialize(buffer);
  data += normVectIndex.serialize(buffer);
  for (const auto &v : normVectArray) {
    data += v.serialize(buffer);
  }

  return data;
}

inline bool BlockIndexKey::operator<(const BlockIndexKey &other) const noexcept {
  if (chrom1 != other.chrom1) {
    return chrom1 < other.chrom1;
  }
  if (chrom2 != other.chrom2) {
    return chrom2 < other.chrom2;
  }
  return resolution < other.resolution;
}

inline HiCFileWriter::HiCFileWriter(HiCHeader header, std::int32_t compression_lvl,
                                    std::size_t buffer_size)
    : _header(std::make_shared<HiCHeader>(std::move(header))),
      _fs(std::make_shared<filestream::FileStream>(
          filestream::FileStream::create(std::string{url()}))),
      _compressor(libdeflate_alloc_compressor(compression_lvl)),
      _compression_buffer(buffer_size, '\0') {
  if (!chromosomes().at(0).is_all()) {
    throw std::runtime_error("reference should contain chromosome ALL");
  }

  _bin_tables.reserve(resolutions().size());
  for (const auto &resolution : resolutions()) {
    _bin_tables.emplace(resolution, BinTableFixed(chromosomes(), resolution));
  }
}

inline std::string_view HiCFileWriter::url() const noexcept {
  assert(_header);
  return _header->url;
}

inline const Reference &HiCFileWriter::chromosomes() const noexcept {
  assert(_header);
  return _header->chromosomes;
}

inline const BinTable &HiCFileWriter::bins(std::uint32_t resolution) const {
  return _bin_tables.at(resolution);
}

inline const std::vector<std::uint32_t> HiCFileWriter::resolutions() const noexcept {
  assert(_header);
  return _header->resolutions;
}

template <typename PixelIt, typename>
inline void HiCFileWriter::append_pixels(std::uint32_t resolution, PixelIt first_pixel,
                                         PixelIt last_pixel) {
  const auto &bin_table = bins(resolution);
  std::for_each(first_pixel, last_pixel, [&](const ThinPixel<float> &p) {
    Pixel<float> pp(bin_table, p);
    const auto key = std::make_pair(pp.coords.bin1.chrom(), pp.coords.bin2.chrom());
    auto it = _pixel_tank.find(key);
    if (it != _pixel_tank.end()) {
      it->second.emplace(p);
    } else {
      _pixel_tank.emplace(key, ChromPixelTank{p});
    }
  });
}

inline void HiCFileWriter::write_pixels() {
  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);
      for (const auto &res : resolutions()) {
        write_pixels(chrom1, chrom2, res);
      }
    }
  }
}

inline void HiCFileWriter::write_pixels(const Chromosome &chrom1, const Chromosome &chrom2,
                                        std::uint32_t resolution) {
  phmap::btree_map<std::uint64_t, std::vector<ThinPixel<float>>> blocks;

  if (!_pixel_tank.contains({chrom1, chrom2})) {
    return;
  }

  {
    const auto &bin_table = bins(resolution);
    assert(!_pixel_tank.empty());
    const auto &pixels = _pixel_tank.at({chrom1, chrom2});
    const auto num_bins = compute_num_bins(chrom1.id(), chrom2.id(), bin_table.bin_size());
    const auto num_columns =
        compute_block_column_count(num_bins, bin_table.bin_size(),
                                   chrom1 == chrom2 ? DEFAULT_INTRA_CUTOFF : DEFAULT_INTER_CUTOFF);
    const auto num_rows = num_bins / num_columns + 1;

    const BlockMapperIntra mapper_intra(num_rows, num_columns);
    const BlockMapperInter mapper_inter(num_rows, num_columns);

    for (const auto &p : pixels) {
      Pixel<float> pp(bin_table, p);

      const auto bin1_id = pp.coords.bin1.rel_id();
      const auto bin2_id = pp.coords.bin2.rel_id();

      const auto block_id =
          pp.coords.is_intra() ? mapper_intra(bin1_id, bin2_id) : mapper_inter(bin1_id, bin2_id);
      auto it = blocks.find(block_id);
      if (it != blocks.end()) {
        it->second.emplace_back(p);
      } else {
        blocks.emplace(block_id, std::vector<ThinPixel<float>>{p});
      }
    }
  }

  for (auto &[block_id, pixels] : blocks) {
    assert(!pixels.empty());
    const auto bin1 = bins(resolution).at(pixels.front().bin1_id);
    const auto bin2 = bins(resolution).at(pixels.front().bin2_id);
    const auto bin1_offset = bin1.rel_id();
    const auto bin2_offset = bin2.rel_id();
    write_interaction_block(block_id, chrom1, chrom2, resolution, pixels, bin1_offset, bin2_offset);
  }
}

inline void HiCFileWriter::write_header() {
  assert(_fs->tellg() == 0);

  assert(_header->version == 9);
  assert(!_header->chromosomes.empty());
  assert(!_header->resolutions.empty());

  _fs->write("HIC\0", 4);
  _fs->write(_header->version);
  _fs->write(std::int64_t(-1));  // masterIndexOffset
  _fs->write(_header->genomeID.c_str(), _header->genomeID.size() + 1);
  _fs->write(_header->nviPosition);
  _fs->write(_header->nviLength);

  // Write attributes
  const auto nAttributes = static_cast<std::int32_t>(_header->attributes.size());
  _fs->write(nAttributes);
  for (const auto &[k, v] : _header->attributes) {
    _fs->write(k.c_str(), k.size() + 1);
    _fs->write(v.c_str(), v.size() + 1);
  }

  // Write chromosomes
  auto numChromosomes = static_cast<std::uint32_t>(_header->chromosomes.size());
  _fs->write(numChromosomes);

  for (const Chromosome &c : _header->chromosomes) {
    const auto name = std::string{c.name()};
    _fs->write(name.c_str(), name.size() + 1);
    _fs->write<std::int64_t>(c.size());
  }

  // write resolutions
  _fs->write(static_cast<std::int32_t>(_header->resolutions.size()));
  const std::vector<std::int32_t> resolutions(_header->resolutions.begin(),
                                              _header->resolutions.end());
  _fs->write(resolutions);
}

inline void HiCFileWriter::write_master_index_offset(std::int64_t master_index_offset) {
  const auto foffset = _fs->tellp();
  const auto offset = sizeof("HIC") + sizeof(_header->version);
  _fs->seekp(offset);
  _fs->write(master_index_offset);
  _fs->seekp(static_cast<std::int64_t>(foffset));
}

inline std::size_t HiCFileWriter::write_matrix_metadata(std::uint32_t chrom1_id,
                                                        std::uint32_t chrom2_id) {
  assert(_header);

  const auto offset = _fs->tellp();

  MatrixMetadata m;

  m.chr1Idx = static_cast<std::int32_t>(chrom1_id);
  m.chr2Idx = static_cast<std::int32_t>(chrom2_id);

  // chrom_id is 0 when chromosome is ALL
  m.nResolutions = chrom1_id == 0 ? std::int32_t(1) : std::int32_t(resolutions().size());

  _fs->write(m.serialize(_bbuffer));

  return offset;
}

inline auto HiCFileWriter::write_resolutions_metadata(std::uint32_t chrom1_id,
                                                      std::uint32_t chrom2_id,
                                                      const std::string &unit) {
  struct Result {
    std::size_t file_offset{};
    std::size_t matrix_metadata_bytes{};
  };

  Result res{};

  // chrom_id is 0 when chromosome is ALL
  const auto nResolutions = chrom1_id == 0 ? std::size_t(1) : std::size_t(resolutions().size());

  const auto offset1 = _fs->tellp();
  for (std::size_t i = 0; i < nResolutions; ++i) {
    const auto bin_size = resolutions()[i];
    const auto num_bins = compute_num_bins(chrom1_id, chrom2_id, bin_size);
    const auto num_columns = compute_block_column_count(
        num_bins, bin_size, chrom1_id == chrom2_id ? DEFAULT_INTRA_CUTOFF : DEFAULT_INTER_CUTOFF);
    const auto num_rows = num_bins / num_columns + 1;

    MatrixResolutionMetadata m{};
    m.unit = unit;
    m.resIdx = static_cast<std::int32_t>(i);
    m.sumCounts = 0;
    m.occupiedCellCount = 0;  // not used
    m.percent5 = 0;           // not used
    m.percent95 = 0;          // not used
    m.binSize = static_cast<std::int32_t>(bin_size);
    m.blockSize = static_cast<std::int32_t>(num_rows);
    m.blockColumnCount = static_cast<std::int32_t>(num_columns);

    const BlockIndexKey key{chromosomes().at(chrom1_id), chromosomes().at(chrom2_id), bin_size};

    auto it = _block_index.find(key);
    m.blockCount = static_cast<std::int32_t>(it != _block_index.end() ? it->second.size() : 0);

    // write resolution metadata
    _fs->write(m.serialize(_bbuffer));

    // write block index

    if (it != _block_index.end()) {
      for (const auto &blk : it->second) {
        _fs->write(blk.serialize(_bbuffer));
      }
    }
  }

  const auto offset2 = _fs->tellp();

  res.file_offset = offset1;
  res.matrix_metadata_bytes = offset2 - offset1;

  return res;
}

inline auto HiCFileWriter::write_body_metadata(const Chromosome &chrom1, const Chromosome &chrom2,
                                               const std::string &unit) {
  struct Result {
    std::size_t matrix_metadata_offset;
    std::size_t matrix_metadata_bytes;
  };

  Result offsets{};
  offsets.matrix_metadata_offset = write_matrix_metadata(chrom1.id(), chrom2.id());
  write_resolutions_metadata(chrom1.id(), chrom2.id(), unit);

  offsets.matrix_metadata_bytes = _fs->tellp() - offsets.matrix_metadata_offset;

  for (const auto &res : resolutions()) {
    _block_index.erase(BlockIndexKey{chrom1, chrom2, res});
  }
  return offsets;
}

inline std::streamoff HiCFileWriter::write_interaction_block(
    std::uint64_t block_id, const Chromosome &chrom1, const Chromosome &chrom2,
    std::uint32_t resolution, const std::vector<ThinPixel<float>> &pixels,
    std::size_t bin_column_offset, std::size_t bin_row_offset) {
  const auto offset = _fs->tellp();

  // TODO support dense layout
  // TODO support representation using shorts

  MatrixInteractionBlock m(bins(resolution), pixels, bin_row_offset);

  m.nRecords = static_cast<std::int32_t>(pixels.size());
  m.binColumnOffset = static_cast<std::int32_t>(bin_column_offset);
  m.binRowOffset = static_cast<std::int32_t>(bin_row_offset);
  m.useFloatContact = true;
  m.useIntXPos = true;
  m.useIntYPos = true;
  m.matrixRepresentation = 1;

  _fs->write(m.serialize(_bbuffer, *_compressor, _compression_buffer));

  MatrixBlockMetadata mm{static_cast<std::int32_t>(block_id), static_cast<std::int64_t>(offset),
                         static_cast<std::int32_t>(_fs->tellp() - offset)};

  const BlockIndexKey key{chrom1, chrom2, resolution};
  auto idx = _block_index.find(key);
  if (idx != _block_index.end()) {
    idx->second.emplace(std::move(mm));
  } else {
    _block_index.emplace(key, phmap::btree_set<MatrixBlockMetadata>{std::move(mm)});
  }
  return static_cast<std::streamoff>(offset);
}

inline std::streamoff HiCFileWriter::write_footers() {
  const auto offset = _fs->tellp();
  // TODO write bytes v5;
  const std::int64_t nBytesV5 = 0;  // TODO
  const std::int32_t nEntries = static_cast<std::int32_t>(_footers.size());
  _fs->write(nBytesV5);
  _fs->write(nEntries);

  for (const auto &f : _footers) {
    _fs->write(f.masterIndex.serialize(_bbuffer));
  }

  for (const auto &f : _footers) {
    _fs->write(f.expectedValues.serialize(_bbuffer));
  }

  for (const auto &f : _footers) {
    _fs->write(f.normExpectedValues.serialize(_bbuffer));
  }

  for (const auto &f : _footers) {
    _fs->write(f.normVectIndex.serialize(_bbuffer));
  }

  for (const auto &f : _footers) {
    for (const auto &v : f.normVectArray) {
      _fs->write(v.serialize(_bbuffer));
    }
  }

  return static_cast<std::streamoff>(offset);
}

inline void HiCFileWriter::add_footer(const Chromosome &chrom1, const Chromosome &chrom2,
                                      std::size_t file_offset, std::size_t matrix_metadata_bytes) {
  FooterV5 f;
  assert(file_offset < _fs->size());
  assert(file_offset + matrix_metadata_bytes <= _fs->size());

  f.masterIndex.key = fmt::format(FMT_STRING("{}_{}"), chrom1.id(), chrom2.id());
  f.masterIndex.position = static_cast<std::int64_t>(file_offset);
  f.masterIndex.size = static_cast<std::int32_t>(matrix_metadata_bytes);

  // TODO populate other fields
  _footers.emplace_back(std::move(f));
}

inline void HiCFileWriter::finalize() {
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, std::size_t> file_offsets{};
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, std::size_t> matrix_size_bytes{};

  std::size_t master_index_offset = std::numeric_limits<std::size_t>::max();

  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);
      const auto offsets = write_body_metadata(chrom1, chrom2);
      const auto key = std::make_pair(chrom1, chrom2);

      master_index_offset = std::min(offsets.matrix_metadata_offset, master_index_offset);
      file_offsets[key] = offsets.matrix_metadata_offset;
      matrix_size_bytes[key] = offsets.matrix_metadata_bytes;
    }
  }

  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);
      const auto key = std::make_pair(chrom1, chrom2);
      add_footer(chrom1, chrom2, file_offsets[key], matrix_size_bytes[key]);
      if (chrom2.is_all()) {
        break;
      }
    }
  }

  const auto footer_offset = write_footers();

  write_master_index_offset(footer_offset);
}

inline std::size_t HiCFileWriter::compute_block_column_count(std::size_t num_bins,
                                                             std::uint32_t bin_size,
                                                             std::uint32_t cutoff,
                                                             std::size_t block_capacity) {
  auto num_columns = num_bins / block_capacity + 1;
  if (bin_size < cutoff) {
    const auto genome_size = num_bins * bin_size;
    num_columns = genome_size / (block_capacity * cutoff);
  }

  const auto max_sqrt =
      static_cast<std::size_t>(std::sqrt(std::numeric_limits<std::int32_t>::max()));
  return std::min(num_columns, max_sqrt - 1);
}

inline std::size_t HiCFileWriter::compute_num_bins(std::uint32_t chrom1_id, std::uint32_t chrom2_id,
                                                   std::size_t bin_size) {
  const auto max_size =
      std::max(chromosomes().at(chrom1_id).size(), chromosomes().at(chrom2_id).size());
  return (max_size + bin_size - 1) / bin_size;
}

inline HiCFileWriter::BlockMapperInter::BlockMapperInter(std::uint64_t block_bin_count,
                                                         std::uint64_t block_column_count)
    : _block_bin_count(block_bin_count), _block_column_count(block_column_count) {
  assert(_block_bin_count != 0);
  assert(_block_column_count != 0);
}

inline std::uint64_t HiCFileWriter::BlockMapperInter::block_column_count() const {
  return _block_column_count;
}

inline std::uint64_t HiCFileWriter::BlockMapperInter::block_bin_count() const {
  return _block_bin_count;
}

inline std::uint64_t HiCFileWriter::BlockMapperInter::operator()(std::uint64_t bin1_id,
                                                                 std::uint64_t bin2_id) const {
  const auto i = bin1_id / block_bin_count();
  const auto j = bin2_id / block_bin_count();

  return (block_column_count() * j) + i;
}

inline HiCFileWriter::BlockMapperIntra::BlockMapperIntra(std::uint64_t block_bin_count,
                                                         std::uint64_t block_column_count,
                                                         std::int64_t base_depth)
    : _inter_mapper(block_bin_count, block_column_count), _base(init_base(base_depth)) {}

inline std::uint64_t HiCFileWriter::BlockMapperIntra::block_column_count() const {
  return _inter_mapper.block_column_count();
}

inline std::uint64_t HiCFileWriter::BlockMapperIntra::block_bin_count() const {
  return _inter_mapper.block_bin_count();
}

inline bool HiCFileWriter::BlockMapperIntra::use_inter_mapper() const noexcept {
  return _base == 0;
}

inline std::uint64_t HiCFileWriter::BlockMapperIntra::operator()(std::uint64_t bin1_id,
                                                                 std::uint64_t bin2_id) const {
  if (use_inter_mapper()) {
    return _inter_mapper(bin1_id, bin2_id);
  }
  const auto delta = bin1_id > bin2_id ? bin1_id - bin2_id : bin2_id - bin1_id;
  const auto n =
      static_cast<double>(delta) / std::sqrt(2.0) / static_cast<double>(block_bin_count());

  const auto depth = static_cast<std::uint64_t>(std::log(1.0 + n) / _base);
  const auto position_along_diagonal = (bin1_id + bin2_id) / 2 / block_bin_count();

  return depth * block_column_count() + position_along_diagonal;
}

inline double HiCFileWriter::BlockMapperIntra::init_base(std::int64_t base_depth) noexcept {
  if (base_depth > 1) {
    return std::log(static_cast<double>(base_depth));
  }
  if (base_depth < 0) {
    return static_cast<double>(-base_depth);
  }
  return std::log(2.0);
}

}  // namespace hictk::hic::internal
