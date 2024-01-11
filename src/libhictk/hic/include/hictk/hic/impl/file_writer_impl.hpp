// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/compile.h>
#include <fmt/format.h>
#include <libdeflate.h>
#include <parallel_hashmap/phmap.h>
#include <spdlog/spdlog.h>
#include <zstd.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <ios>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/filestream.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/reference.hpp"
#include "hictk/transformers/coarsen.hpp"

namespace hictk::hic::internal {

inline bool BlockIndexKey::operator<(const BlockIndexKey &other) const noexcept {
  if (chrom1 != other.chrom1) {
    return chrom1 < other.chrom1;
  }
  if (chrom2 != other.chrom2) {
    return chrom2 < other.chrom2;
  }
  return resolution < other.resolution;
}

inline bool MatrixBodyMetadataTank::Key::operator==(
    const MatrixBodyMetadataTank::Key &other) const noexcept {
  return chrom1 == other.chrom1 && chrom2 == other.chrom2;
}

inline bool MatrixBodyMetadataTank::contains(const Chromosome &chrom1,
                                             const Chromosome &chrom2) const noexcept {
  return _tank.contains(Key{chrom1, chrom2});
}

inline auto MatrixBodyMetadataTank::at(const Chromosome &chrom1, const Chromosome &chrom2) const
    -> const MatrixBodyMetadata & {
  return _tank.at({chrom1, chrom2});
}

inline HiCSectionOffsets MatrixBodyMetadataTank::offset(const Chromosome &chrom1,
                                                        const Chromosome &chrom2) const {
  return _offsets.at(Key{chrom1, chrom2});
}

inline void MatrixBodyMetadataTank::insert(const Chromosome &chrom1, const Chromosome &chrom2,
                                           MatrixMetadata matrix_metadata,
                                           MatrixResolutionMetadata matrix_resolution_metadata) {
  auto match = _tank.find(Key{chrom1, chrom2});
  if (match != _tank.end()) {
    match->second.matrixMetadata = std::move(matrix_metadata);
    match->second.resolutionMetadata.push_back(std::move(matrix_resolution_metadata));
  } else {
    _tank.emplace(
        Key{chrom1, chrom2},
        MatrixBodyMetadata{std::move(matrix_metadata), std::vector<MatrixResolutionMetadata>{
                                                           std::move(matrix_resolution_metadata)}});
  }
}

inline void MatrixBodyMetadataTank::update_offsets(const Chromosome &chrom1,
                                                   const Chromosome &chrom2,
                                                   std::streamoff position, std::size_t size) {
  _offsets.emplace(Key{chrom1, chrom2}, HiCSectionOffsets{position, size});
}

inline void MatrixBodyMetadataTank::remove(const Chromosome &chrom1, const Chromosome &chrom2) {
  const Key k{chrom1, chrom2};
  _tank.erase(k);
  _offsets.erase(k);
}

inline auto MatrixBodyMetadataTank::operator()() const noexcept
    -> const phmap::flat_hash_map<Key, MatrixBodyMetadata> & {
  return _tank;
}

inline std::string MatrixBodyMetadataTank::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  for (const auto &[_, metadata] : _tank) {
    std::ignore = metadata.serialize(buffer, false);
  }

  return buffer.get();
}

inline HiCFileWriter::HiCFileWriter(HiCHeader header, const std::filesystem::path &tmpdir,
                                    std::int32_t compression_lvl, std::size_t buffer_size)
    : _header(std::make_shared<const HiCHeader>(std::move(header))),
      _fs(std::make_shared<filestream::FileStream>(filestream::FileStream::create(_header->url))),
      _tmpdir(tmpdir.empty() ? nullptr
                             : std::make_unique<const hictk::internal::TmpDir>(
                                   tmpdir / (_header->url + ".tmp/"))),
      _bin_tables(init_bin_tables(chromosomes(), resolutions())),
      _block_mappers(init_interaction_block_mappers((*_tmpdir)(), _bin_tables, 3)),
      _compressor(libdeflate_alloc_compressor(compression_lvl)),
      _compression_buffer(buffer_size, '\0') {}

inline std::string_view HiCFileWriter::url() const noexcept {
  assert(_header);
  return _header->url;
}

inline const Reference &HiCFileWriter::chromosomes() const noexcept {
  assert(_header);
  return _header->chromosomes;
}

inline const BinTable &HiCFileWriter::bins(std::uint32_t resolution) const {
  return *_bin_tables.at(resolution);
}

inline const std::vector<std::uint32_t> &HiCFileWriter::resolutions() const noexcept {
  assert(_header);
  return _header->resolutions;
}

inline void HiCFileWriter::write_header() {
  assert(_fs->tellp() == 0);

  assert(_header);
  assert(_header->version == 9);
  assert(!chromosomes().empty());

  const auto offset1 = _fs->tellp();

  _fs->write("HIC\0", 4);
  _fs->write(_header->version);
  _fs->write(std::int64_t(-1));  // masterIndexOffset
  _fs->write(_header->genomeID.c_str(), _header->genomeID.size() + 1);
  _fs->write(std::int64_t(_fs->tellp() + sizeof(_header->normVectorIndexPosition)));
  _fs->write(std::int64_t(0));

  // Write attributes
  const auto nAttributes = static_cast<std::int32_t>(_header->attributes.size());
  _fs->write(nAttributes);
  for (const auto &[k, v] : _header->attributes) {
    _fs->write(k.c_str(), k.size() + 1);
    _fs->write(v.c_str(), v.size() + 1);
  }

  // Write chromosomes
  auto numChromosomes = static_cast<std::uint32_t>(chromosomes().size());
  _fs->write(numChromosomes);

  for (const Chromosome &c : chromosomes()) {
    const auto name = std::string{c.name()};
    _fs->write(name.c_str(), name.size() + 1);
    _fs->write<std::int64_t>(c.size());
  }

  // write resolutions
  _fs->write(static_cast<std::int32_t>(_header->resolutions.size()));
  const std::vector<std::int32_t> resolutions(_header->resolutions.begin(),
                                              _header->resolutions.end());
  _fs->write(resolutions);

  // write fragments: TODO
  const std::int32_t nFragResolutions = 0;
  _fs->write(nFragResolutions);

  const auto offset2 = _fs->tellp();

  _header_offsets.position = static_cast<std::streamoff>(offset1);
  _header_offsets.size = offset2 - offset1;
}

inline void HiCFileWriter::write_footer_offset() {
  const auto foffset = _fs->tellp();
  const auto offset = sizeof("HIC") + sizeof(_header->version);
  _fs->seekp(offset);
  _fs->write(conditional_static_cast<std::int64_t>(_footer_offsets.position));
  _fs->seekp(static_cast<std::int64_t>(foffset));
}

inline void HiCFileWriter::write_norm_vector_index(std::streamoff position, std::size_t length) {
  const auto foffset = _fs->tellp();
  const auto offset =
      static_cast<std::int64_t>(sizeof("HIC") + sizeof(_header->version) +
                                sizeof(_header->masterIndexOffset) + _header->genomeID.size() + 1);
  _fs->seekp(offset);
  _fs->write(static_cast<std::int64_t>(position));
  _fs->write(static_cast<std::int64_t>(length));
  _fs->seekp(static_cast<std::int64_t>(foffset));
}

template <typename PixelIt, typename>
inline void HiCFileWriter::add_pixels(PixelIt first_pixel, PixelIt last_pixel) {
  add_pixels(resolutions().front(), first_pixel, last_pixel);
}

inline auto HiCFileWriter::write_pixels(const Chromosome &chrom1, const Chromosome &chrom2)
    -> HiCSectionOffsets {
  const auto pos = _fs->tellp();
  auto block_section = write_pixels(chrom1, chrom2, resolutions().front());
  add_body_metadata(resolutions().front(), chrom1, chrom2);
  write_body_metadata();
  add_footer(chrom1, chrom2);
  _footer_offsets = write_footers();

  finalize();

  for (std::size_t i = 1; i < resolutions().size(); ++i) {
    auto base_resolution = resolutions().front();
    const auto res = resolutions()[i];

    for (std::size_t j = 0; j < i; ++j) {
      if (resolutions()[j] % res == 0) {
        base_resolution = resolutions()[j];
      }
    }

    const File f(std::string{url()}, base_resolution);
    const auto sel = f.fetch();
    const auto factor = res / base_resolution;
    const transformers::CoarsenPixels coarsener(
        sel.begin<float>(), sel.end<float>(),
        std::make_shared<const BinTable>(bins(base_resolution)), factor);

    auto &mapper = _block_mappers.at(res);
    mapper.append_pixels(coarsener.begin(), coarsener.end());

    const auto end_of_block_section =
        block_section.position + static_cast<std::streamoff>(block_section.size);
    _fs->seekp(end_of_block_section);
    block_section = write_pixels(chrom1, chrom2, res);
    for (std::size_t j = 0; j < i; ++j) {
      add_body_metadata(resolutions()[j], chrom1, chrom2);
    }
    _body_metadata_offsets = write_body_metadata();
    add_footer(chrom1, chrom2);
    _footer_offsets = write_footers();
    finalize();
  }
  return {static_cast<std::streamoff>(pos), _fs->tellp() - pos};
}

inline auto HiCFileWriter::write_body_metadata() -> HiCSectionOffsets {
  const auto pos = _fs->tellp();
  for (const auto &[chroms, metadata] : _matrix_metadata()) {
    const auto pos1 = _fs->tellp();
    _fs->write(metadata.serialize(_bbuffer));
    const auto pos2 = _fs->tellp();
    _matrix_metadata.update_offsets(chroms.chrom1, chroms.chrom2, static_cast<std::streamoff>(pos1),
                                    pos2 - pos1);
  }
  const auto size = _fs->tellp() - static_cast<std::size_t>(pos);

  return {static_cast<std::streamoff>(pos), size};
}

inline void HiCFileWriter::add_body_metadata(std::uint32_t resolution, const Chromosome &chrom1,
                                             const Chromosome &chrom2, const std::string &unit) {
  const auto sum_counts = _block_mappers.at(resolution).pixel_sum(chrom1, chrom2);
  if (sum_counts == 0) {
    return;
  }

  auto metadata = _matrix_metadata.contains(chrom1, chrom2) ? _matrix_metadata.at(chrom1, chrom2)
                                                            : MatrixBodyMetadata{};

  auto &mm = metadata.matrixMetadata;
  MatrixResolutionMetadata mrm{};

  const auto num_bins = compute_num_bins(chrom1, chrom2, resolution);
  const auto num_columns = compute_block_column_count(chrom1, chrom2, resolution);
  const auto num_rows = num_bins / num_columns + 1;

  mrm.unit = unit;
  mrm.resIdx = static_cast<std::int32_t>(std::distance(
      resolutions().begin(), std::find(resolutions().begin(), resolutions().end(), resolution)));
  mrm.sumCounts = sum_counts;
  mrm.occupiedCellCount = 0;  // not used
  mrm.percent5 = 0;           // not used
  mrm.percent95 = 0;          // not used
  mrm.binSize = static_cast<std::int32_t>(resolution);
  mrm.blockSize = static_cast<std::int32_t>(num_rows);
  mrm.blockColumnCount = static_cast<std::int32_t>(num_columns);

  metadata.resolutionMetadata.emplace_back(std::move(mrm));
  metadata.matrixMetadata.chr1Idx = static_cast<std::int32_t>(chrom1.id());
  metadata.matrixMetadata.chr2Idx = static_cast<std::int32_t>(chrom2.id());
  metadata.matrixMetadata.nResolutions =
      static_cast<std::int32_t>(metadata.resolutionMetadata.size());

  _matrix_metadata.insert(chrom1, chrom2, mm, mrm);
}

inline auto HiCFileWriter::write_footers() -> HiCSectionOffsets {
  const auto offset1 = _fs->tellp();
  std::int64_t nBytesV5 = -1;
  const auto nEntries = static_cast<std::int32_t>(_footers.size());
  _fs->write(nBytesV5);
  _fs->write(nEntries);

  for (const auto &[_, footer] : _footers) {
    _fs->write(footer.serialize(_bbuffer));
  }

  const auto offset2 = _fs->tellp();
  nBytesV5 = static_cast<std::int64_t>(offset2 - offset1);

  _fs->seekp(static_cast<std::streamoff>(offset1));
  _fs->write(nBytesV5);
  _fs->seekp(static_cast<std::streamoff>(offset2));

  const auto normVectorIndexPosition = _fs->tellp();
  _fs->write(std::int32_t(0));  // no nNormExpectedValueVectors
  _fs->write(std::int32_t(0));  // no NormVectors
  const auto normVectorIndexLength = _fs->tellp() - normVectorIndexPosition;

  write_norm_vector_index(static_cast<std::streamoff>(normVectorIndexPosition),
                          normVectorIndexLength);

  return {static_cast<std::streamoff>(offset1), offset2 - offset1};
}

inline void HiCFileWriter::add_footer(const Chromosome &chrom1, const Chromosome &chrom2) {
  FooterV5 footer{};
  if (!_matrix_metadata.contains(chrom1, chrom2)) {
    return;
  }
  footer.masterIndex.key = fmt::format(FMT_STRING("{}_{}"), chrom1.id(), chrom2.id());

  const auto offset = _matrix_metadata.offset(chrom1, chrom2);
  footer.masterIndex.position = static_cast<std::int64_t>(offset.position);
  footer.masterIndex.size = static_cast<std::int32_t>(offset.size);

  _footers.emplace(std::make_pair(chrom1, chrom2), std::move(footer));
}

inline void HiCFileWriter::finalize() {
  write_footer_offset();
  _fs->flush();
}

inline auto HiCFileWriter::init_bin_tables(const Reference &chromosomes,
                                           const std::vector<std::uint32_t> &resolutions)
    -> BinTables {
  BinTables bin_tables(resolutions.size());
  for (const auto &res : resolutions) {
    bin_tables.emplace(res, std::make_shared<const BinTable>(chromosomes, res));
  }

  return bin_tables;
}

inline auto HiCFileWriter::init_interaction_block_mappers(const std::filesystem::path &root_folder,
                                                          const BinTables &bin_tables,
                                                          int compression_lvl) -> BlockMappers {
  BlockMappers mappers(bin_tables.size());
  for (const auto &[res, bin_table] : bin_tables) {
    const auto path = fmt::format(FMT_STRING("{}/{}.bin"), root_folder.string(), res);
    mappers.emplace(res, HiCInteractionToBlockMapper{path, bin_table, compression_lvl});
  }

  return mappers;
}

template <typename PixelIt, typename>
inline void HiCFileWriter::add_pixels(std::uint32_t resolution, PixelIt first_pixel,
                                      PixelIt last_pixel) {
  auto &mapper = _block_mappers.at(resolution);
  mapper.append_pixels(first_pixel, last_pixel);
}

inline auto HiCFileWriter::write_pixels(const Chromosome &chrom1, const Chromosome &chrom2,
                                        std::uint32_t resolution) -> HiCSectionOffsets {
  const auto offset = _fs->tellp();

  auto &mapper = _block_mappers.at(resolution);
  mapper.finalize();
  const auto &block_ids = mapper.chromosome_index().at(std::make_pair(chrom1, chrom2));
  for (const auto &bid : block_ids) {
    const auto pixels = mapper.merge_blocks(bid);
    write_interaction_block(bid.bid, chrom1, chrom2, resolution, pixels);
  }

  return {static_cast<std::streamoff>(offset), _fs->tellp() - offset};
}

inline auto HiCFileWriter::write_interaction_block(std::uint64_t block_id, const Chromosome &chrom1,
                                                   const Chromosome &chrom2,
                                                   std::uint32_t resolution,
                                                   const MatrixInteractionBlock<float> &blk)
    -> HiCSectionOffsets {
  const auto offset = _fs->tellp();

  _fs->write(blk.serialize(_bbuffer, *_compressor, _compression_buffer));

  MatrixBlockMetadata mm{static_cast<std::int32_t>(block_id), static_cast<std::int64_t>(offset),
                         static_cast<std::int32_t>(_fs->tellp() - offset)};

  const BlockIndexKey key{chrom1, chrom2, resolution};
  auto idx = _block_index.find(key);
  if (idx != _block_index.end()) {
    idx->second.emplace(std::move(mm));
  } else {
    _block_index.emplace(key, phmap::btree_set<MatrixBlockMetadata>{std::move(mm)});
  }
  return {static_cast<std::streamoff>(offset), _fs->tellp() - offset};
}

inline std::size_t HiCFileWriter::compute_num_bins(const Chromosome &chrom1,
                                                   const Chromosome &chrom2,
                                                   std::uint32_t resolution) {
  return HiCInteractionToBlockMapper::compute_num_bins(chrom1, chrom2, resolution);
}

inline std::size_t HiCFileWriter::compute_block_column_count(const Chromosome &chrom1,
                                                             const Chromosome &chrom2,
                                                             std::uint32_t resolution) {
  return HiCInteractionToBlockMapper::compute_block_column_count(
      chrom1, chrom2, resolution,
      chrom1 == chrom2 ? HiCInteractionToBlockMapper::DEFAULT_INTRA_CUTOFF
                       : HiCInteractionToBlockMapper::DEFAULT_INTER_CUTOFF);
}

/*
inline ChromChromHiCFileWriter::ChromChromHiCFileWriter(const Chromosome &chrom1,
                                                        const Chromosome &chrom2, HiCHeader header,
                                                        const std::filesystem::path &tmpdir,
                                                        std::int32_t compression_lvl,
                                                        std::size_t buffer_size)
    : _header(init_header(chrom1, chrom2, std::move(header))),
      _fs(std::make_shared<filestream::FileStream>(
          filestream::FileStream::create(std::string{url()}))),
      _bin_table(
          std::make_shared<const BinTable>(_header->chromosomes, _header->resolutions.front())),
      _compressor(libdeflate_alloc_compressor(compression_lvl)),
      _compression_buffer(buffer_size, '\0'),
      _pixel_tank(_bin_table),
      _tmpdir(tmpdir.empty() ? nullptr
                             : std::make_unique<const hictk::internal::TmpDir>(
                                   tmpdir / (_header->url + ".tmp/"))) {
  _chrom1 = chromosomes().at(0);
  _chrom2 = chromosomes().size() == 2 ? chromosomes().at(1) : _chrom1;

  if (_header->resolutions.size() != 1) {
    throw std::runtime_error(
        "ChromChromHiCFileWriter should be constructed with an header containing a single "
        "resolution");
  }
}

inline std::string_view ChromChromHiCFileWriter::url() const noexcept {
  assert(_header);
  return _header->url;
}

inline const Reference &ChromChromHiCFileWriter::chromosomes() const noexcept {
  assert(_header);
  return _header->chromosomes;
}

inline const BinTable &ChromChromHiCFileWriter::bins() const { return *_bin_table; }

inline std::uint32_t ChromChromHiCFileWriter::resolution() const noexcept {
  assert(_header);
  return _header->resolutions.front();
}

template <typename PixelIt, typename>
inline void ChromChromHiCFileWriter::write_pixels(PixelIt first_pixel, PixelIt last_pixel,
                                                  bool update_expected_values) {
  SPDLOG_DEBUG(FMT_STRING("[{}::{}] writing pixels..."),
               std::filesystem::path{_header->url}.filename(), resolution());
  _pixel_tank.add_pixels(first_pixel, last_pixel, update_expected_values);
  _pixel_tank.finalize();
  _body_offset.position =
      _header_offset.position + static_cast<std::streamoff>(_header_offset.size);
  _body_metadata_offset.position = _body_offset.position;

  struct PixelBlock {
    std::vector<ThinPixel<float>> pixels;
    std::uint32_t bin1_id_offset;
    std::uint32_t bin2_id_offset;
  };

  if (!_pixel_tank.contains(_chrom1, _chrom2)) {
    _body_offset.size = 0;
    _body_metadata_offset.size = 0;
    return;
  }

  phmap::btree_map<std::uint64_t, PixelBlock> blocks;

  {
    const auto &bin_table = bins();
    const auto &pixels = _pixel_tank.pixels(_chrom1, _chrom2);
    const auto num_bins = compute_num_bins(_chrom1.id(), _chrom2.id(), bin_table.bin_size());
    const auto num_columns = compute_block_column_count(
        num_bins, bin_table.bin_size(),
        _chrom1 == _chrom2 ? DEFAULT_INTRA_CUTOFF : DEFAULT_INTER_CUTOFF);
    const auto num_rows = num_bins / num_columns + 1;

    SPDLOG_DEBUG(FMT_STRING("[{}::{}] matrix {}:{}: blockBinCount={}"),
                 std::filesystem::path{_header->url}.filename(), resolution(), _chrom1.name(),
                 _chrom2.name(), num_rows);
    SPDLOG_DEBUG(FMT_STRING("[{}::{}] matrix {}:{}: blockColumnCount={}"),
                 std::filesystem::path{_header->url}.filename(), resolution(), _chrom1.name(),
                 _chrom2.name(), num_columns);

    const BlockMapperIntra mapper_intra(num_rows, num_columns);
    const BlockMapperInter mapper_inter(num_rows, num_columns);

    for (const auto &p : pixels) {
      Pixel<float> pp(bin_table, p);

      assert(pp.coords.bin1.chrom() == _chrom1);
      assert(pp.coords.bin2.chrom() == _chrom2);

      const auto &bin1 = pp.coords.bin1;
      const auto &bin2 = pp.coords.bin2;

      const auto bin1_id = bin1.rel_id();
      const auto bin2_id = bin2.rel_id();

      const auto block_id =
          pp.coords.is_intra() ? mapper_intra(bin1_id, bin2_id) : mapper_inter(bin1_id, bin2_id);
      auto it = blocks.find(block_id);
      if (it != blocks.end()) {
        it->second.pixels.emplace_back(p);
        it->second.bin1_id_offset = std::min(it->second.bin1_id_offset, bin1_id);
        it->second.bin2_id_offset = std::min(it->second.bin2_id_offset, bin2_id);
      } else {
        blocks.emplace(block_id, PixelBlock{{p}, bin1_id, bin2_id});
      }
    }
  }

  for (auto &[block_id, blk] : blocks) {
    SPDLOG_DEBUG(FMT_STRING("[{}::{}] writing block {}..."),
                 std::filesystem::path{_header->url}.filename(), resolution(), block_id);
    write_interaction_block(block_id, blk.pixels, blk.bin1_id_offset, blk.bin2_id_offset);
  }

  _body_metadata_offset.position = static_cast<std::streamoff>(_fs->tellp());
  write_body_metadata();
  _body_offset.size = _fs->tellp() - static_cast<std::size_t>(_body_offset.position);
  _body_metadata_offset.size =
      _fs->tellp() - static_cast<std::size_t>(_body_metadata_offset.position);
}

inline void ChromChromHiCFileWriter::write_header() {
  assert(_fs->tellp() == 0);

  assert(_header);
  assert(_header->version == 9);
  assert(!chromosomes().empty());

  const auto offset1 = _fs->tellp();

  _fs->write("HIC\0", 4);
  _fs->write(_header->version);
  _fs->write(std::int64_t(-1));  // masterIndexOffset
  _fs->write(_header->genomeID.c_str(), _header->genomeID.size() + 1);
  _fs->write(std::int64_t(_fs->tellp() + sizeof(_header->normVectorIndexPosition)));
  _fs->write(std::int64_t(0));

  // Write attributes
  const auto nAttributes = static_cast<std::int32_t>(_header->attributes.size());
  _fs->write(nAttributes);
  for (const auto &[k, v] : _header->attributes) {
    _fs->write(k.c_str(), k.size() + 1);
    _fs->write(v.c_str(), v.size() + 1);
  }

  // Write chromosomes
  auto numChromosomes = static_cast<std::uint32_t>(chromosomes().size());
  _fs->write(numChromosomes);

  for (const Chromosome &c : chromosomes()) {
    const auto name = std::string{c.name()};
    _fs->write(name.c_str(), name.size() + 1);
    _fs->write<std::int64_t>(c.size());
  }

  // write resolutions
  _fs->write(static_cast<std::int32_t>(_header->resolutions.size()));
  const std::vector<std::int32_t> resolutions(_header->resolutions.begin(),
                                              _header->resolutions.end());
  _fs->write(resolutions);

  // write fragments: TODO
  const std::int32_t nFragResolutions = 0;
  _fs->write(nFragResolutions);

  const auto offset2 = _fs->tellp();

  _header_offset.position = static_cast<std::streamoff>(offset1);
  _header_offset.size = offset2 - offset1;
}

inline void ChromChromHiCFileWriter::write_footer_offset(std::int64_t master_index_offset) {
  const auto foffset = _fs->tellp();
  const auto offset = sizeof("HIC") + sizeof(_header->version);
  _fs->seekp(offset);
  _fs->write(master_index_offset);
  _fs->seekp(static_cast<std::int64_t>(foffset));
}

inline void ChromChromHiCFileWriter::write_norm_vector_index(std::streamoff position,
                                                             std::size_t length) {
  const auto foffset = _fs->tellp();
  const auto offset =
      static_cast<std::int64_t>(sizeof("HIC") + sizeof(_header->version) +
                                sizeof(_header->masterIndexOffset) + _header->genomeID.size() + 1);
  _fs->seekp(offset);
  _fs->write(static_cast<std::int64_t>(position));
  _fs->write(static_cast<std::int64_t>(length));
  _fs->seekp(static_cast<std::int64_t>(foffset));
}

inline auto ChromChromHiCFileWriter::write_matrix_metadata() -> HiCSectionOffsets {
  SPDLOG_DEBUG(FMT_STRING("[{}::{}] writing matrix {}:{} metadata..."),
               std::filesystem::path{_header->url}.filename(), resolution(), _chrom1.name(),
               _chrom2.name());

  const auto offset = _fs->tellp();

  _matrix_metadata.chr1Idx = static_cast<std::int32_t>(_chrom1.id());
  _matrix_metadata.chr2Idx = static_cast<std::int32_t>(_chrom2.id());

  _matrix_metadata.nResolutions = 1;

  _fs->write(_matrix_metadata.serialize(_bbuffer));

  return {static_cast<std::streamoff>(offset), _fs->tellp() - offset};
}

inline auto ChromChromHiCFileWriter::write_resolution_metadata(float sum_counts,
                                                               const std::string &unit)
    -> HiCSectionOffsets {
  SPDLOG_DEBUG(FMT_STRING("[{}::{}] writing matrix resolution {}:{} metadata..."),
               std::filesystem::path{_header->url}.filename(), resolution(), _chrom1.name(),
               _chrom2.name());
  const auto offset1 = _fs->tellp();
  const auto num_bins = compute_num_bins(_chrom1.id(), _chrom2.id(), resolution());
  const auto num_columns = compute_block_column_count(
      num_bins, resolution(), _chrom1 == _chrom2 ? DEFAULT_INTRA_CUTOFF : DEFAULT_INTER_CUTOFF);
  const auto num_rows = num_bins / num_columns + 1;

  MatrixResolutionMetadata m{};
  m.unit = unit;
  m.resIdx = 0;
  m.sumCounts = sum_counts;
  m.occupiedCellCount = 0;  // not used
  m.percent5 = 0;           // not used
  m.percent95 = 0;          // not used
  m.binSize = static_cast<std::int32_t>(resolution());
  m.blockSize = static_cast<std::int32_t>(num_rows);
  m.blockColumnCount = static_cast<std::int32_t>(num_columns);

  const BlockIndexKey key{_chrom1, _chrom2, resolution()};

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

  const auto offset2 = _fs->tellp();

  return {static_cast<std::streamoff>(offset1), offset2 - offset1};
}

inline auto ChromChromHiCFileWriter::write_body_metadata(const std::string &unit)
    -> HiCSectionOffsets {
  if (!_pixel_tank.contains(_chrom1, _chrom2)) {
    return {static_cast<std::streamoff>(_fs->tellp()), 0};
  }
  const auto tot_counts = _pixel_tank.matrix_counts(_chrom1, _chrom2);

  const auto pos = _fs->tellp();
  write_matrix_metadata();
  write_resolution_metadata(tot_counts, unit);
  const auto size = _fs->tellp() - static_cast<std::size_t>(pos);

  return {static_cast<std::streamoff>(pos), size};
}

inline auto ChromChromHiCFileWriter::write_interaction_block(
    std::uint64_t block_id, const std::vector<ThinPixel<float>> &pixels,
    std::size_t bin_column_offset, std::size_t bin_row_offset) -> HiCSectionOffsets {
  const auto offset = _fs->tellp();

  // TODO support dense layout
  // TODO support representation using shorts

  MatrixInteractionBlock m(bins(), pixels, bin_row_offset);

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

  const BlockIndexKey key{_chrom1, _chrom2, resolution()};
  auto idx = _block_index.find(key);
  if (idx != _block_index.end()) {
    idx->second.emplace(std::move(mm));
  } else {
    _block_index.emplace(key, phmap::btree_set<MatrixBlockMetadata>{std::move(mm)});
  }
  return {static_cast<std::streamoff>(offset), _fs->tellp() - offset};
}

inline auto ChromChromHiCFileWriter::write_footer() -> HiCSectionOffsets {
  _footer.masterIndex.key = fmt::format(FMT_STRING("{}_{}"), _chrom1.id(), _chrom2.id());

  _footer.masterIndex.position = static_cast<std::int64_t>(_body_metadata_offset.position);
  _footer.masterIndex.size = static_cast<std::int32_t>(_body_metadata_offset.size);

  const auto offset = _fs->tellp();
  const std::int64_t nBytesV5 = -1;
  const std::int32_t nEntries = 1;
  _fs->write(nBytesV5);
  _fs->write(nEntries);

  _fs->write(_footer.masterIndex.serialize(_bbuffer));

  return {static_cast<std::streamoff>(offset), _fs->tellp() - offset};
}

inline void ChromChromHiCFileWriter::write_footer_section_size(std::streamoff footer_offset,
                                                               std::uint64_t bytes) {
  const auto offset = _fs->tellp();
  _fs->seekp(static_cast<std::int64_t>(footer_offset));
  _fs->write(static_cast<std::int64_t>(bytes));
  _fs->seekp(static_cast<std::int64_t>(offset));
}

inline void ChromChromHiCFileWriter::write_expected_values(std::string_view unit) {
  ExpectedValues evs{};
  evs.nExpectedValueVectors = 1;

  std::vector<std::uint32_t> chrom_ids{};
  std::vector<double> scale_factors{};
  const auto &ev = _pixel_tank.expected_values();
  assert(!ev.weights().empty());
  for (const auto &[chrom, factor] : ev.scaling_factors()) {
    chrom_ids.push_back(chrom.id());
    scale_factors.push_back(factor);
  }
  evs.expectedValues.emplace_back(unit, resolution(), ev.weights(), chrom_ids, scale_factors);

  _fs->write(evs.serialize(_bbuffer));
}

inline void ChromChromHiCFileWriter::finalize() {
  const auto [footer_pos, footer_size] = write_footer();

  write_expected_values("BP");

  write_footer_section_size(footer_pos, _fs->tellp() - static_cast<std::size_t>(footer_pos));

  const auto normVectorIndexPosition = _fs->tellp();
  _fs->write(std::int32_t(0));  // no nNormExpectedValueVectors
  _fs->write(std::int32_t(0));  // no NormVectors
  const auto normVectorIndexLength = _fs->tellp() - normVectorIndexPosition;

  write_footer_offset(footer_pos);
  write_norm_vector_index(static_cast<std::streamoff>(normVectorIndexPosition),
                          normVectorIndexLength);
}
 */

}  // namespace hictk::hic::internal
