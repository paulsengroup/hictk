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

template <typename N>
inline PixelTank<N>::PixelTank(std::shared_ptr<const BinTable> bins)
    : _bin_table(std::move(bins)), _expected_values(_bin_table) {}

template <typename N>
inline bool PixelTank<N>::contains(const Chromosome &chrom1,
                                   const Chromosome &chrom2) const noexcept {
  return _pixel_tank.contains(std::make_pair(chrom1, chrom2));
}

template <typename N>
inline void PixelTank<N>::add_pixel(const ThinPixel<N> &p, bool update_expected_values) {
  assert(_bin_table);
  add_pixel(Pixel<N>{*_bin_table, p}, update_expected_values);
}

template <typename N>
inline void PixelTank<N>::add_pixel(const Pixel<N> &p, bool update_expected_values) {
  const auto key = std::make_pair(p.coords.bin1.chrom(), p.coords.bin2.chrom());
  auto it = _pixel_tank.find(key);
  if (it != _pixel_tank.end()) {
    it->second.emplace_back(p.to_thin());
    assert(_matrix_tot_counts.contains(key));
    _matrix_tot_counts.at(key) += 1;
  } else {
    _pixel_tank.emplace(key, ChromPixelTank{p.to_thin()});
    _matrix_tot_counts.emplace(key, 1.0F);
  }

  if (update_expected_values) {
    _expected_values.add(p);
  }
}

template <typename N>
template <typename PixelIt, typename>
inline void PixelTank<N>::add_pixels(PixelIt first_pixel, PixelIt last_pixel,
                                     bool update_expected_values) {
  using PixelT = decltype(*first_pixel);
  using PixelN = decltype(std::declval<PixelT>().count);
  static_assert(std::is_same_v<PixelN, N>);

  std::for_each(first_pixel, last_pixel,
                [&](const auto &p) { add_pixel(p, update_expected_values); });
}

template <typename N>
inline void PixelTank<N>::finalize() {
  SPDLOG_DEBUG(FMT_STRING("[{}] finalize called on PixelTank"), _bin_table->bin_size());
  for (auto &[k, v] : _pixel_tank) {
    SPDLOG_DEBUG(FMT_STRING("[{}] sorting {} pixels for PixelTank[{}:{}]"), _bin_table->bin_size(),
                 v.size(), k.first.name(), k.second.name());
    std::sort(v.begin(), v.end());  // TODO is this needed?
  }
  _expected_values.compute_density();
}

template <typename N>
inline auto PixelTank<N>::pixels() const noexcept
    -> const phmap::flat_hash_map<ChromosomePair, ChromPixelTank> & {
  return _pixel_tank;
}

template <typename N>
inline auto PixelTank<N>::matrix_counts() const noexcept
    -> const phmap::flat_hash_map<ChromosomePair, float> & {
  return _matrix_tot_counts;
}

template <typename N>
inline auto PixelTank<N>::expected_values() const noexcept {
  return _expected_values;
}

template <typename N>
inline auto PixelTank<N>::pixels(const Chromosome &chrom1, const Chromosome &chrom2) const
    -> const ChromPixelTank & {
  return _pixel_tank.at(std::make_pair(chrom1, chrom2));
}

template <typename N>
inline auto PixelTank<N>::matrix_counts(const Chromosome &chrom1, const Chromosome &chrom2) const
    -> const float & {
  return _matrix_tot_counts.at(std::make_pair(chrom1, chrom2));
}

inline bool MetadataOffsetTank::Key::operator<(
    const MetadataOffsetTank::Key &other) const noexcept {
  if (resolution != other.resolution) {
    return resolution < other.resolution;
  }
  if (chrom1 != other.chrom1) {
    return chrom1 < other.chrom1;
  }
  return chrom2 < other.chrom2;
}

inline bool MetadataOffsetTank::contains(const Chromosome &chrom1, const Chromosome &chrom2,
                                         std::uint32_t resolution) const noexcept {
  return _tank.contains({chrom1, chrom2, resolution});
}

inline auto MetadataOffsetTank::at(const Chromosome &chrom1, const Chromosome &chrom2,
                                   std::uint32_t resolution) const -> const Value & {
  return _tank.at({chrom1, chrom2, resolution});
}

inline void MetadataOffsetTank::insert(const Chromosome &chrom1, const Chromosome &chrom2,
                                       std::uint32_t resolution, std::size_t offset,
                                       std::size_t size) {
  SPDLOG_DEBUG(FMT_STRING("[{}] body offsets for {}:{}: {} + {}"), resolution, chrom1.name(),
               chrom2.name(), offset, size);
  assert(!contains(chrom1, chrom2, resolution));
  _tank.emplace(Key{chrom1, chrom2, resolution}, Value{static_cast<std::streamoff>(offset), size});
}

inline auto MetadataOffsetTank::operator()() const noexcept
    -> const phmap::btree_map<Key, Value> & {
  return _tank;
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

inline std::size_t ChromChromHiCFileWriter::compute_block_column_count(std::size_t num_bins,
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
  return std::clamp(num_columns, std::size_t(1), max_sqrt - 1);
}

inline std::shared_ptr<const HiCHeader> ChromChromHiCFileWriter::init_header(
    const Chromosome &chrom1, const Chromosome &chrom2, HiCHeader &&header) {
  header.chromosomes = chrom1 != chrom2 ? Reference{chrom1, chrom2} : Reference{chrom1};

  return std::make_shared<const HiCHeader>(std::move(header));
}

inline std::size_t ChromChromHiCFileWriter::compute_num_bins(std::uint32_t chrom1_id,
                                                             std::uint32_t chrom2_id,
                                                             std::size_t bin_size) {
  const auto max_size =
      std::max(chromosomes().at(chrom1_id).size(), chromosomes().at(chrom2_id).size());
  return (max_size + bin_size - 1) / bin_size;
}
 */

}  // namespace hictk::hic::internal
