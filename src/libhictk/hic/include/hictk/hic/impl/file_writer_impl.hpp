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

  BinaryBuffer tmp_buffer;
  for (const auto &b : blocksMetadata) {
    if (b.blockSizeBytes != 0) {
      buffer.write(b.serialize(tmp_buffer));
    }
  }

  return buffer.get();
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
  std::string data{};

  for (const auto &idx : masterIndex) {
    buffer.clear();
    data += idx.serialize(buffer);
  }
  buffer.clear();
  data += expectedValues.serialize(buffer);
  buffer.clear();
  data += normExpectedValues.serialize(buffer);
  buffer.clear();
  data += normVectIndex.serialize(buffer);
  for (const auto &v : normVectArray) {
    buffer.clear();
    data += v.serialize(buffer);
  }

  const auto total_size = static_cast<std::int64_t>(sizeof(nBytesV5) + sizeof(nEntries)) +
                          static_cast<std::int64_t>(data.size() * sizeof(char));

  buffer.clear();
  buffer.write(total_size);
  buffer.write(nEntries);

  return buffer.get() + data;
}

constexpr std::int64_t FooterV5::master_index_offset() const noexcept { return sizeof(nBytesV5); }

inline void FooterV5::add_footer(const HiCFooter &footer, std::int64_t matrix_offset,
                                 std::int32_t matrix_bytes) {
  ++nEntries;
  masterIndex.emplace_back(
      MasterIndex{fmt::format(FMT_STRING("{}_{}\0"), footer.chrom1().id(), footer.chrom2().id()),
                  matrix_offset, matrix_bytes}

  );
  // TODO populate other fields
}

inline HiCFileWriter::HiCFileWriter(HiCHeader header)
    : _header(std::make_shared<HiCHeader>(std::move(header))),
      _fs(std::make_shared<filestream::FileStream>(
          filestream::FileStream::create(std::string{url()}))) {
  assert(_header->chromosomes.contains("ALL"));
}

inline std::string_view HiCFileWriter::url() const noexcept {
  assert(_header);
  return _header->url;
}

inline const Reference &HiCFileWriter::chromosomes() const noexcept {
  assert(_header);
  return _header->chromosomes;
}

inline const std::vector<std::uint32_t> HiCFileWriter::resolutions() const noexcept {
  assert(_header);
  return _header->resolutions;
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

  // const auto genome_size = _header->chromosomes.chrom_size_prefix_sum().back();

  // _fs->write("ALL", 4);
  // _fs->write(static_cast<std::int64_t>(genome_size));

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

inline std::pair<std::streamoff, std::size_t> HiCFileWriter::write_body_metadata(
    std::uint32_t chrom1_id, std::uint32_t chrom2_id, const std::string &unit) {
  const auto offset = _fs->tellp();
  write_matrix_metadata(chrom1_id, chrom2_id);
  write_resolution_metadata(chrom1_id, chrom2_id, unit);

  _block_index.clear();

  return std::make_pair(static_cast<std::streamoff>(offset), _fs->tellp() - offset);
}

inline void HiCFileWriter::write_matrix_metadata(std::uint32_t chrom1_id, std::uint32_t chrom2_id) {
  assert(_header);

  // chrom_id is 0 when chromosome is ALL
  const auto nResolutions = chrom1_id == 0 ? std::int32_t(1) : std::int32_t(resolutions().size());

  _fs->write(static_cast<std::int32_t>(chrom1_id));
  _fs->write(static_cast<std::int32_t>(chrom2_id));
  _fs->write(nResolutions);
}

inline void HiCFileWriter::write_resolution_metadata(std::uint32_t chrom1_id,
                                                     std::uint32_t chrom2_id,
                                                     const std::string &unit) {
  // chrom_id is 0 when chromosome is ALL
  const auto nResolutions = chrom1_id == 0 ? std::size_t(1) : std::size_t(resolutions().size());

  for (std::size_t i = 0; i < nResolutions; ++i) {
    const auto bin_size = resolutions()[i];
    const auto num_bins = compute_num_bins(chrom1_id, chrom2_id, bin_size);
    const auto num_columns = compute_block_column_count(
        num_bins, bin_size, chrom1_id == chrom2_id ? DEFAULT_INTRA_CUTOFF : DEFAULT_INTER_CUTOFF);
    const auto num_rows = num_bins / (num_columns + 1);

    // write resolution metadata
    const auto resIdx = static_cast<std::int32_t>(i);
    const float sumCount = 0;
    const float occupiedCellCount = 0;
    const float percent5 = 0;
    const float percent95 = 0;
    const auto binSize = static_cast<std::int32_t>(bin_size);
    const std::int32_t blockSize = static_cast<std::int32_t>(num_rows);
    const auto blockColumnCount = static_cast<std::int32_t>(num_columns);
    const auto blockCount = static_cast<std::int32_t>(_block_index.size());

    _fs->write(unit.c_str(), unit.size() + 1);
    _fs->write(resIdx);
    _fs->write(sumCount);
    _fs->write(occupiedCellCount);
    _fs->write(percent5);
    _fs->write(percent95);
    _fs->write(binSize);
    _fs->write(blockSize);
    _fs->write(blockColumnCount);
    _fs->write(blockCount);

    // write block index
    for (const auto &b : _block_index) {
      _fs->write(b.serialize(_bbuffer));
    }
  }
}

inline std::streamoff HiCFileWriter::write_interaction_block(const InteractionBlock &blk,
                                                             std::size_t bin_column_offset,
                                                             std::size_t bin_row_offset) {
  const auto offset = _fs->tellp();
  _bbuffer.clear();

  // TODO support dense layout
  // TODO support representation using shorts

  const auto nRecords = static_cast<std::int32_t>(blk.size());
  const auto binColumnOffset = static_cast<std::int32_t>(bin_column_offset);
  const auto binRowOffset = static_cast<std::int32_t>(bin_row_offset);
  const auto useFloatContact = true;
  const auto useIntXPos = true;
  const auto useIntYPos = true;
  const std::uint8_t matrixRepresentation = 1;

  _bbuffer.write(nRecords);
  _bbuffer.write(binColumnOffset);
  _bbuffer.write(binRowOffset);
  _bbuffer.write(useFloatContact);
  _bbuffer.write(useIntXPos);
  _bbuffer.write(useIntYPos);
  _bbuffer.write(matrixRepresentation);

  const auto interactions = group_interactions_by_column(blk, bin_row_offset);
  const auto rowCount = static_cast<std::int32_t>(interactions.size());  // TODO support short
  _bbuffer.write(rowCount);

  for (const auto &[row, pixels] : interactions) {
    const auto rowNumber = static_cast<std::int32_t>(row);              // TODO support short
    const auto recordCount = static_cast<std::int32_t>(pixels.size());  // TODO support short
    _bbuffer.write(rowNumber);
    _bbuffer.write(recordCount);

    assert(std::is_sorted(pixels.begin(), pixels.end()));
    for (const auto &p : pixels) {
      const auto binColumn = static_cast<std::int32_t>(p.bin1_id - bin_row_offset);
      const auto value = p.count;
      _bbuffer.write(binColumn);
      _bbuffer.write(value);
    }
  }

  // TODO compress properly
  auto *compressor = libdeflate_alloc_compressor(9);
  std::string out(1'000'000, '\0');
  const auto compressed_size = libdeflate_zlib_compress(
      compressor, _bbuffer.get().data(), _bbuffer.get().size(), out.data(), out.size());
  libdeflate_free_compressor(compressor);
  out.resize(compressed_size);
  _fs->write(out);

  MatrixBlockMetadata m{static_cast<std::int32_t>(blk.id()), static_cast<std::int64_t>(offset),
                        static_cast<std::int32_t>(_bbuffer.get().size())};

  assert(!_block_index.contains(m));
  _block_index.emplace(std::move(m));

  return static_cast<std::streamoff>(offset);
}

inline std::streamoff HiCFileWriter::write_footer(const std::vector<HiCFooter> &footers,
                                                  const std::vector<std::int64_t> &matrix_offsets,
                                                  const std::vector<std::int32_t> &matrix_bytes) {
  assert(footers.size() == matrix_offsets.size());
  assert(footers.size() == matrix_bytes.size());

  const auto offset = _fs->tellp();

  FooterV5 footer{};

  for (std::size_t i = 0; i < footers.size(); ++i) {
    footer.add_footer(footers[i], matrix_offsets[i], matrix_bytes[i]);
  }

  _fs->write(footer.serialize(_bbuffer));

  return static_cast<std::streamoff>(offset);
}

inline phmap::btree_map<std::int32_t, std::vector<ThinPixel<float>>>
HiCFileWriter::group_interactions_by_column(const InteractionBlock &blk,
                                            std::size_t bin_row_offset) {
  phmap::btree_map<std::int32_t, std::vector<ThinPixel<float>>> buffer;

  for (const auto &p : blk) {
    const auto col = static_cast<std::int32_t>(p.bin2_id - bin_row_offset);
    auto it = buffer.find(col);
    if (it != buffer.end()) {
      it->second.push_back(p);
    } else {
      buffer.emplace(col, std::vector<ThinPixel<float>>{p});
    }
  }

  for (auto &[_, v] : buffer) {
    std::sort(v.begin(), v.end());
  }
  return buffer;
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
  return max_size / (bin_size + 1);
}

inline HiCFileWriter::BlockMapperInter::BlockMapperInter(std::uint64_t block_bin_count,
                                                         std::uint64_t block_column_count)
    : _block_bin_count(block_bin_count), _block_column_count(block_column_count) {}

inline std::uint64_t HiCFileWriter::BlockMapperInter::block_column_count() {
  return _block_column_count;
}

inline std::uint64_t HiCFileWriter::BlockMapperInter::block_bin_count() { return _block_bin_count; }

inline std::uint64_t HiCFileWriter::BlockMapperInter::operator()(std::uint64_t bin1_id,
                                                                 std::uint64_t bin2_id) {
  const auto i = bin1_id / block_bin_count();
  const auto j = bin2_id / block_bin_count();

  return (block_column_count() * j) + i;
}

inline HiCFileWriter::BlockMapperIntra::BlockMapperIntra(std::uint64_t block_bin_count,
                                                         std::uint64_t block_column_count,
                                                         std::int64_t base_depth)
    : _inter_mapper(block_bin_count, block_column_count), _base(init_base(base_depth)) {}

inline std::uint64_t HiCFileWriter::BlockMapperIntra::block_column_count() {
  return _inter_mapper.block_column_count();
}

inline std::uint64_t HiCFileWriter::BlockMapperIntra::block_bin_count() {
  return _inter_mapper.block_bin_count();
}

inline bool HiCFileWriter::BlockMapperIntra::use_inter_mapper() const noexcept {
  return _base == 0;
}

inline std::uint64_t HiCFileWriter::BlockMapperIntra::operator()(std::uint64_t bin1_id,
                                                                 std::uint64_t bin2_id) {
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
