// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#if __has_include(<blockingconcurrentqueue.h>)
#include <blockingconcurrentqueue.h>
#else
#include <concurrentqueue/blockingconcurrentqueue.h>
#endif
#include <fmt/compile.h>
#include <fmt/format.h>
#include <libdeflate.h>
#include <parallel_hashmap/phmap.h>
#if __has_include(<readerwriterqueue.h>)
#include <readerwriterqueue.h>
#else
#include <readerwriterqueue/readerwriterqueue.h>
#endif
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
#include <queue>
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

template <typename I1, typename I2>
inline HiCSectionOffsets::HiCSectionOffsets(I1 start_, I2 size_)
    : _position(conditional_static_cast<std::streamoff>(start_)),
      _size(conditional_static_cast<std::size_t>(size_)) {
  static_assert(std::is_integral_v<I1>);
  static_assert(std::is_integral_v<I2>);
}

inline std::streamoff HiCSectionOffsets::start() const noexcept { return _position; }

inline std::streamoff HiCSectionOffsets::end() const noexcept {
  return _position + static_cast<std::streamoff>(size());
}

inline std::size_t HiCSectionOffsets::size() const noexcept { return _size; }

inline std::size_t &HiCSectionOffsets::size() noexcept { return _size; }

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
    match->second.resolutionMetadata.emplace(std::move(matrix_resolution_metadata));
  } else {
    _tank.emplace(
        Key{chrom1, chrom2},
        MatrixBodyMetadata{std::move(matrix_metadata), phmap::btree_set<MatrixResolutionMetadata>{
                                                           std::move(matrix_resolution_metadata)}});
  }
}

inline void MatrixBodyMetadataTank::update_offsets(const Chromosome &chrom1,
                                                   const Chromosome &chrom2,
                                                   std::streamoff position, std::size_t size) {
  auto [it, inserted] = _offsets.emplace(Key{chrom1, chrom2}, HiCSectionOffsets{position, size});
  if (!inserted) {
    it->second = HiCSectionOffsets{position, size};
  }
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

inline HiCFileWriter::HiCFileWriter(HiCHeader header, std::size_t n_threads,
                                    const std::filesystem::path &tmpdir,
                                    std::int32_t compression_lvl, std::size_t buffer_size)
    : _header(init_header(std::move(header))),
      _fs(std::make_shared<filestream::FileStream>(filestream::FileStream::create(_header->url))),
      _tmpdir(tmpdir.empty() ? nullptr
                             : std::make_unique<const hictk::internal::TmpDir>(
                                   tmpdir / (_header->url + ".tmp/"))),
      _bin_tables(init_bin_tables(chromosomes(), resolutions())),
      _block_mappers(init_interaction_block_mappers((*_tmpdir)(), _bin_tables, 3)),
      _compression_lvl(compression_lvl),
      _compressor(libdeflate_alloc_compressor(compression_lvl)),
      _compression_buffer(buffer_size, '\0'),
      _tpool(init_tpool(n_threads)) {}

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

inline void HiCFileWriter::serialize() {
  write_header();
  write_pixels();
  compute_and_write_expected_values();
  finalize();
}

inline void HiCFileWriter::write_header() {
  assert(_fs->tellp() == 0);

  assert(_header);
  assert(_header->version == 9);
  assert(!chromosomes().empty());

  const auto offset1 = _fs->tellp();

  SPDLOG_DEBUG(FMT_STRING("writing header at offset {}"), offset1);

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

  _header_section = {offset1, offset2 - offset1};
  _data_block_section = {offset2, 0};
  _body_metadata_section = {offset2, 0};
  _footer_section = {offset2, 0};
}

inline void HiCFileWriter::write_footer_size() {
  SPDLOG_DEBUG(FMT_STRING("updating footer size to {}"), _footer_section.size());
  const auto nBytesV5 = static_cast<std::int64_t>(_footer_section.size());
  _fs->seekp(_footer_section.start());
  _fs->write(nBytesV5);
}

inline void HiCFileWriter::write_footer_offset() {
  SPDLOG_DEBUG(FMT_STRING("updating footer offset to {}"), _footer_section.start());
  const auto offset = sizeof("HIC") + sizeof(_header->version);
  _fs->seekp(offset);
  _fs->write(conditional_static_cast<std::int64_t>(_footer_section.start()));
}

inline void HiCFileWriter::write_norm_vector_index() {
  const auto offset =
      static_cast<std::int64_t>(sizeof("HIC") + sizeof(_header->version) +
                                sizeof(_header->footerPosition) + _header->genomeID.size() + 1);
  const auto normVectorIndexPosition =
      conditional_static_cast<std::int64_t>(_expected_values_norm_section.start());
  const auto normVectorIndexLength = static_cast<std::int64_t>(
      _expected_values_norm_section.size() + _norm_vectors_section.size());

  SPDLOG_DEBUG(FMT_STRING("writing normVectorIndex {}:{} at offset {}..."), normVectorIndexPosition,
               normVectorIndexLength, offset);
  _fs->seekp(offset);
  _fs->write(normVectorIndexPosition);
  _fs->write(normVectorIndexLength);
}

template <typename PixelIt, typename>
inline void HiCFileWriter::add_pixels(PixelIt first_pixel, PixelIt last_pixel) {
  add_pixels(resolutions().front(), first_pixel, last_pixel);
}

inline void HiCFileWriter::write_pixels() {
  SPDLOG_DEBUG(FMT_STRING("begin writing interaction blocks..."));
  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);
      if (chrom2.is_all()) {
        break;
      }
      write_pixels(chrom1, chrom2);
    }
  }
  write_all_matrix();
}

inline void HiCFileWriter::write_all_matrix(std::uint32_t target_resolution) {
  auto base_resolution = resolutions().front();
  const auto factor = target_resolution / base_resolution;
  target_resolution = factor * base_resolution;

  for (const auto &res : resolutions()) {
    if (res > target_resolution) {
      break;
    }

    if (target_resolution % res == 0) {
      base_resolution = res;
    }
  }

  SPDLOG_DEBUG(FMT_STRING("writing pixels for {}:{} matrix..."), chromosomes().at(0).name(),
               chromosomes().at(0).name());

  // Dummy bin table that always map to chromosome #1
  BinTable bin_table_(
      Reference{Chromosome{0, "__ALL__", std::numeric_limits<std::uint32_t>::max()}}, 1);

  File f(std::string{url()}, base_resolution);
  auto sel = f.fetch();

  MatrixInteractionBlock<float> blk{};
  if (base_resolution == target_resolution) {
    std::for_each(sel.begin<float>(), sel.end<float>(),
                  [&](const ThinPixel<float> &p) { blk.emplace_back(Pixel(bin_table_, p)); });
  } else {
    transformers::CoarsenPixels coarsener(sel.begin<float>(), sel.end<float>(),
                                          std::make_shared<const BinTable>(bins(base_resolution)),
                                          factor);
    std::for_each(coarsener.begin(), coarsener.end(),
                  [&](const ThinPixel<float> &p) { blk.emplace_back(Pixel(bin_table_, p)); });
  }
  const auto chrom = chromosomes().at(0);
  assert(chrom.is_all());
  write_interaction_block(0, chrom, chrom, target_resolution, blk);

  add_body_metadata(target_resolution, chrom, chrom);
  write_body_metadata();
  add_footer(chrom, chrom);
  write_footers();

  finalize();
}

inline auto HiCFileWriter::write_pixels(const Chromosome &chrom1, const Chromosome &chrom2)
    -> HiCSectionOffsets {
  write_pixels(chrom1, chrom2, resolutions().front());
  add_body_metadata(resolutions().front(), chrom1, chrom2);
  write_body_metadata();
  add_footer(chrom1, chrom2);
  write_footers();

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
    const auto sel = f.fetch(chrom1.name(), chrom2.name());
    const auto factor = res / base_resolution;
    const transformers::CoarsenPixels coarsener(
        sel.begin<float>(), sel.end<float>(),
        std::make_shared<const BinTable>(bins(base_resolution)), factor);

    auto &mapper = _block_mappers.at(res);
    mapper.append_pixels(coarsener.begin(), coarsener.end(), _tpool);

    write_pixels(chrom1, chrom2, res);
    for (std::size_t j = 0; j <= i; ++j) {
      add_body_metadata(resolutions()[j], chrom1, chrom2);
    }
    write_body_metadata();
    add_footer(chrom1, chrom2);
    write_footers();
    finalize();
    mapper.clear();
  }
  return {_data_block_section.start(),
          _fs->tellp() - static_cast<std::size_t>(_data_block_section.start())};
}

inline void HiCFileWriter::write_body_metadata() {
  const auto pos = _data_block_section.end();
  _fs->seekp(pos);
  for (const auto &[chroms, metadata] : _matrix_metadata()) {
    const auto pos1 = _fs->tellp();
    SPDLOG_DEBUG(FMT_STRING("writing MatrixBodyMetadata for {}:{} ({} resolutions) at offset {}"),
                 chroms.chrom1.name(), chroms.chrom2.name(), metadata.resolutionMetadata.size(),
                 pos1);
    _fs->write(metadata.serialize(_bbuffer));
    const auto pos2 = _fs->tellp();
    SPDLOG_DEBUG(FMT_STRING("updating MatrixBodyMetadata offset and size for {}:{} ({} "
                            "resolutions) to {} and {}"),
                 chroms.chrom1.name(), chroms.chrom2.name(), metadata.resolutionMetadata.size(),
                 pos1, pos2 - pos1);
    _matrix_metadata.update_offsets(chroms.chrom1, chroms.chrom2, static_cast<std::streamoff>(pos1),
                                    pos2 - pos1);
  }

  const auto size = _fs->tellp() - static_cast<std::size_t>(pos);
  _body_metadata_section = {pos, size};
}

inline void HiCFileWriter::add_body_metadata(std::uint32_t resolution, const Chromosome &chrom1,
                                             const Chromosome &chrom2, const std::string &unit) {
  SPDLOG_DEBUG(FMT_STRING("adding MatrixBodyMetadata for {}:{} at {} {}"), chrom1.name(),
               chrom2.name(), resolution, unit);
  const auto sum_counts =
      chrom1.is_all() ? 1.0F : _block_mappers.at(resolution).pixel_sum(chrom1, chrom2);
  if (sum_counts == 0) {
    return;
  }

  auto metadata = _matrix_metadata.contains(chrom1, chrom2) ? _matrix_metadata.at(chrom1, chrom2)
                                                            : MatrixBodyMetadata{};

  auto &mm = metadata.matrixMetadata;
  MatrixResolutionMetadata mrm{};

  const auto num_bins = compute_num_bins(chrom1, chrom2, resolution);
  const auto num_columns =
      chrom1.is_all() ? std::size_t(1) : compute_block_column_count(chrom1, chrom2, resolution);
  const auto num_rows = chrom1.is_all() ? num_bins : num_bins / num_columns + 1;

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

  const auto &blks = _block_index.at(BlockIndexKey{chrom1, chrom2, resolution});
  mrm.set_block_metadata(blks.begin(), blks.end());

  mm.chr1Idx = static_cast<std::int32_t>(chrom1.id());
  mm.chr2Idx = static_cast<std::int32_t>(chrom2.id());
  mm.nResolutions = static_cast<std::int32_t>(metadata.resolutionMetadata.size() + 1);

  _matrix_metadata.insert(chrom1, chrom2, mm, mrm);
}

inline void HiCFileWriter::write_footers() {
  const auto offset1 = _body_metadata_section.end();
  _fs->seekp(offset1);
  SPDLOG_DEBUG(FMT_STRING("initializing footer section at offset {}"), offset1);
  const std::int64_t nBytesV5 = -1;
  const auto nEntries = static_cast<std::int32_t>(_footers.size());
  _fs->write(nBytesV5);
  _fs->write(nEntries);

  for (const auto &[chroms, footer] : _footers) {
    SPDLOG_DEBUG(FMT_STRING("writing FooterV5 for {}:{} at offset {}"), chroms.first.name(),
                 chroms.second.name(), _fs->tellp());
    _fs->write(footer.serialize(_bbuffer));
  }

  write_empty_expected_values();

  _footer_section = {offset1, _fs->tellp() - static_cast<std::size_t>(offset1)};
}

inline void HiCFileWriter::add_footer(const Chromosome &chrom1, const Chromosome &chrom2) {
  if (!_matrix_metadata.contains(chrom1, chrom2)) {
    return;
  }

  FooterV5 footer{};
  footer.masterIndex.key = fmt::format(FMT_STRING("{}_{}"), chrom1.id(), chrom2.id());

  const auto offset = _matrix_metadata.offset(chrom1, chrom2);
  footer.masterIndex.position = conditional_static_cast<std::int64_t>(offset.start());
  footer.masterIndex.size = static_cast<std::int32_t>(offset.size());

  auto [it, inserted] = _footers.emplace(std::make_pair(chrom1, chrom2), footer);
  if (!inserted) {
    it->second = std::move(footer);
  }
}

inline void HiCFileWriter::write_empty_expected_values() {
  ExpectedValues ev{};
  ev.nExpectedValueVectors = 0;

  const auto offset = _fs->tellp();
  SPDLOG_DEBUG(FMT_STRING("writing empty expected values section at offset {}..."), offset);
  _fs->write(ev.serialize(_bbuffer));

  _expected_values_section = {offset, _fs->tellp() - offset};
}

inline void HiCFileWriter::write_empty_normalized_expected_values() {
  const auto offset = _expected_values_section.end();
  SPDLOG_DEBUG(FMT_STRING("writing empty expected values (normalized) section at offset {}..."),
               offset);
  _fs->seekp(offset);
  _fs->write(std::int32_t(0));
  _expected_values_norm_section = {offset, _fs->tellp() - static_cast<std::size_t>(offset)};
}

inline void HiCFileWriter::write_empty_norm_vectors() {
  const auto offset = _expected_values_norm_section.end();
  SPDLOG_DEBUG(FMT_STRING("writing empty normalization vector section at offset {}..."), offset);
  _fs->seekp(offset);
  _fs->write(std::int32_t(0));
  _norm_vectors_section = {offset, _fs->tellp() - static_cast<std::size_t>(offset)};
}

inline void HiCFileWriter::compute_and_write_expected_values() {
  ExpectedValues ev{};
  ev.nExpectedValueVectors = static_cast<std::int32_t>(resolutions().size());

  for (const auto &resolution : resolutions()) {
    SPDLOG_DEBUG(FMT_STRING("computing expected values at resolution {}..."), resolution);
    File f(std::string{url()}, resolution);
    auto sel = f.fetch();

    ExpectedValuesAggregator aggr(_bin_tables.at(resolution));
    std::for_each(sel.begin<float>(), sel.end<float>(), [&](const auto &p) { aggr.add(p); });
    aggr.compute_density();

    std::vector<float> weights(aggr.weights().size());
    std::transform(aggr.weights().begin(), aggr.weights().end(), weights.begin(),
                   [](const auto w) { return static_cast<float>(w); });

    std::vector<std::uint32_t> chrom_ids{};
    std::vector<double> scaling_factors{};
    std::for_each(aggr.scaling_factors().begin(), aggr.scaling_factors().end(),
                  [&](const auto &kv) {
                    chrom_ids.push_back(kv.first.id());
                    scaling_factors.push_back(kv.second);
                  });

    ev.expectedValues.emplace_back(
        ExpectedValuesBlock{"BP", resolution, aggr.weights(), chrom_ids, scaling_factors});
  }

  const auto offset =
      _footer_section.end() - static_cast<std::streamoff>(sizeof(ev.nExpectedValueVectors));
  _fs->seekp(offset);
  _fs->write(ev.serialize(_bbuffer));

  _expected_values_section = {offset, _fs->tellp() - static_cast<std::size_t>(offset)};
  _footer_section.size() += _expected_values_section.size() - sizeof(ev.nExpectedValueVectors);
}

inline void HiCFileWriter::finalize() {
  write_footer_size();
  write_footer_offset();
  write_empty_normalized_expected_values();
  write_empty_norm_vectors();
  write_norm_vector_index();
  _fs->flush();
  _fs->seekp(0, std::ios::end);
}

inline std::shared_ptr<const HiCHeader> HiCFileWriter::init_header(HiCHeader &&header,
                                                                   std::uint32_t all_scale_factor) {
  header.chromosomes = header.chromosomes.add_ALL(all_scale_factor);
  return std::make_shared<const HiCHeader>(std::move(header));
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

inline BS::thread_pool HiCFileWriter::init_tpool(std::size_t n_threads) {
  return {conditional_static_cast<BS::concurrency_t>(n_threads < 2 ? std::size_t(1) : n_threads)};
}

template <typename PixelIt, typename>
inline void HiCFileWriter::add_pixels(std::uint32_t resolution, PixelIt first_pixel,
                                      PixelIt last_pixel) {
  auto &mapper = _block_mappers.at(resolution);
  mapper.append_pixels(first_pixel, last_pixel, _tpool);
}

inline auto HiCFileWriter::write_pixels(const Chromosome &chrom1, const Chromosome &chrom2,
                                        std::uint32_t resolution) -> HiCSectionOffsets {
  SPDLOG_DEBUG(FMT_STRING("writing pixels for {}:{} matrix at {} resolution..."), chrom1.name(),
               chrom2.name(), resolution);

  const auto offset = _data_block_section.end();
  _fs->seekp(offset);

  const auto pixels_written = write_interaction_blocks(chrom1, chrom2, resolution);

  SPDLOG_DEBUG(FMT_STRING("written {} pixels for {}:{} matrix at {} resolution"), pixels_written,
               chrom1.name(), chrom2.name(), resolution);

  _data_block_section.size() += _fs->tellp() - static_cast<std::size_t>(offset);
  return {offset, _fs->tellp() - static_cast<std::size_t>(offset)};
}

inline std::size_t HiCFileWriter::write_interaction_blocks(const Chromosome &chrom1,
                                                           const Chromosome &chrom2,
                                                           std::uint32_t resolution) {
  auto &mapper = _block_mappers.at(resolution);
  mapper.finalize();

  const auto block_ids = mapper.chromosome_index().find(std::make_pair(chrom1, chrom2));
  if (block_ids == mapper.chromosome_index().end()) {
    SPDLOG_DEBUG(FMT_STRING("no pixels to write for {}:{} matrix at {} resolution"), chrom1.name(),
                 chrom2.name(), resolution);
    return 0;
  }

  if (_tpool.get_thread_count() < 3 || block_ids->second.size() == 1) {
    std::size_t pixels_written = 0;
    for (const auto &bid : block_ids->second) {
      auto blk = mapper.merge_blocks(bid);
      pixels_written += static_cast<std::size_t>(blk.nRecords);
      write_interaction_block(bid.bid, chrom1, chrom2, resolution, std::move(blk));
    }

    return pixels_written;
  }

  std::mutex block_id_queue_mtx{};
  std::queue<std::uint64_t> block_id_queue{};
  moodycamel::BlockingConcurrentQueue<std::pair<std::uint64_t, MatrixInteractionBlock<float>>>
      block_queue(_tpool.get_thread_count() * 4);

  std::mutex serialized_block_tank_mtx{};
  phmap::flat_hash_map<std::uint64_t, std::string> serialized_block_tank{_tpool.get_thread_count() *
                                                                         4};
  const auto stop_token = std::numeric_limits<std::uint64_t>::max();
  std::atomic<bool> early_return = false;

  std::mutex mapper_mtx{};

  std::vector<std::future<void>> worker_threads{};
  for (BS::concurrency_t i = 2; i < _tpool.get_thread_count(); ++i) {
    worker_threads.emplace_back(_tpool.submit([&]() {
      compress_blocks_thr(block_queue, serialized_block_tank, serialized_block_tank_mtx,
                          early_return, stop_token);
    }));
  }

  auto producer = _tpool.submit([&]() {
    return merge_blocks_thr(block_ids->second, mapper, block_queue, block_id_queue,
                            block_id_queue_mtx, mapper_mtx, early_return,
                            std::vector<std::uint64_t>(worker_threads.size(), stop_token));
  });

  auto writer = _tpool.submit([&]() {
    write_compressed_blocks_thr(chrom1, chrom2, resolution, block_id_queue, block_id_queue_mtx,
                                serialized_block_tank, serialized_block_tank_mtx, early_return,
                                stop_token);
  });

  const auto pixels_written = producer.get();
  for (auto &worker : worker_threads) {
    worker.get();
  }
  writer.get();

  return pixels_written;
}

inline auto HiCFileWriter::write_interaction_block(std::uint64_t block_id, const Chromosome &chrom1,
                                                   const Chromosome &chrom2,
                                                   std::uint32_t resolution,
                                                   const MatrixInteractionBlock<float> &blk)
    -> HiCSectionOffsets {
  const auto offset = _fs->tellp();

  std::ignore = blk.serialize(_bbuffer, *_compressor, _compression_buffer);
  SPDLOG_DEBUG(FMT_STRING("writing block #{} for {}:{}:{} at {}:{}"), block_id, chrom1.name(),
               chrom2.name(), resolution, offset, _compression_buffer.size());
  _fs->write(_compression_buffer);

  MatrixBlockMetadata mm{static_cast<std::int32_t>(block_id), static_cast<std::int64_t>(offset),
                         static_cast<std::int32_t>(_fs->tellp() - offset)};

  const BlockIndexKey key{chrom1, chrom2, resolution};
  auto idx = _block_index.find(key);
  if (idx != _block_index.end()) {
    idx->second.emplace(std::move(mm));
  } else {
    _block_index.emplace(key, phmap::btree_set<MatrixBlockMetadata>{std::move(mm)});
  }
  return {offset, _fs->tellp() - offset};
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

inline std::size_t HiCFileWriter::merge_blocks_thr(
    const phmap::btree_set<HiCInteractionToBlockMapper::BlockID> &block_ids,
    HiCInteractionToBlockMapper &mapper,
    moodycamel::BlockingConcurrentQueue<std::pair<std::uint64_t, MatrixInteractionBlock<float>>>
        &block_queue,
    std::queue<std::uint64_t> &block_id_queue, std::mutex &block_id_queue_mtx,
    std::mutex &mapper_mtx, std::atomic<bool> &early_return,
    const std::vector<std::uint64_t> &stop_tokens) {
  try {
    std::size_t pixels_written = 0;
    for (const auto &bid : block_ids) {
      if (early_return) {
        return pixels_written;
      }

      // merge blocks corresponding to bid
      auto blk = mapper.merge_blocks(bid, mapper_mtx);
      pixels_written += static_cast<std::size_t>(blk.nRecords);

      // enqueue block
      while (!block_queue.try_enqueue(std::make_pair(bid.bid, blk))) {
        std::this_thread::sleep_for(std::chrono::milliseconds(50));
        if (early_return) {
          return pixels_written;
        }
      }
      std::scoped_lock lck(block_id_queue_mtx);
      block_id_queue.push(bid.bid);
    }

    // signal to consumers that no more blocks will be enqueued
    for (const auto &tok : stop_tokens) {
      block_queue.enqueue(std::make_pair(tok, MatrixInteractionBlock<float>{}));
      std::scoped_lock lck(block_id_queue_mtx);
      block_id_queue.push(tok);
    }

    return pixels_written;
  } catch (...) {
    early_return = true;
    throw;
  }
}

inline void HiCFileWriter::compress_blocks_thr(
    moodycamel::BlockingConcurrentQueue<std::pair<std::uint64_t, MatrixInteractionBlock<float>>>
        &block_queue,
    phmap::flat_hash_map<std::uint64_t, std::string> &serialized_block_tank,
    std::mutex &serialized_block_tank_mtx, std::atomic<bool> &early_return,
    std::uint64_t stop_token) {
  try {
    std::pair<std::uint64_t, MatrixInteractionBlock<float>> buffer{};
    BinaryBuffer bbuffer{};
    std::string compression_buffer(16'000'000, '\0');
    std::unique_ptr<libdeflate_compressor> compressor(
        libdeflate_alloc_compressor(_compression_lvl));

    while (!early_return) {
      // dequeue block
      if (!block_queue.wait_dequeue_timed(buffer, std::chrono::milliseconds(500))) {
        continue;
      }
      const auto &[bid, blk] = buffer;
      if (bid == stop_token) {
        return;
      }

      // compress and serialize block
      std::ignore = blk.serialize(bbuffer, *compressor, compression_buffer);

      // enqueue serialized block
      std::scoped_lock lck(serialized_block_tank_mtx);
      serialized_block_tank.emplace(std::make_pair(bid, compression_buffer));
    }
  } catch (...) {
    early_return = true;
    throw;
  }
}

inline void HiCFileWriter::write_compressed_blocks_thr(
    const Chromosome &chrom1, const Chromosome &chrom2, std::uint32_t resolution,
    std::queue<std::uint64_t> &block_id_queue, std::mutex &block_id_queue_mtx,
    phmap::flat_hash_map<std::uint64_t, std::string> &serialized_block_tank,
    std::mutex &serialized_block_tank_mtx, std::atomic<bool> &early_return,
    std::uint64_t stop_token) {
  try {
    std::string buffer;
    while (!early_return) {
      if (block_id_queue.empty()) {
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
        continue;
      }

      const auto bid = [&]() {
        std::scoped_lock lck(block_id_queue_mtx);
        const auto n = block_id_queue.front();
        block_id_queue.pop();
        return n;
      }();

      if (bid == stop_token) {
        return;
      }
      while (!early_return) {
        {
          std::scoped_lock lck(serialized_block_tank_mtx);
          auto match = serialized_block_tank.find(bid);
          if (match != serialized_block_tank.end()) {
            buffer = match->second;
            serialized_block_tank.erase(match);
            break;
          }
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(10));
      }
      if (early_return) {
        return;
      }

      const auto offset = _fs->tellp();
      SPDLOG_DEBUG(FMT_STRING("writing block #{} for {}:{}:{} at {}:{}"), bid, chrom1.name(),
                   chrom2.name(), resolution, offset, buffer.size());
      _fs->write(buffer);

      MatrixBlockMetadata mm{static_cast<std::int32_t>(bid), static_cast<std::int64_t>(offset),
                             static_cast<std::int32_t>(_fs->tellp() - offset)};
      const BlockIndexKey key{chrom1, chrom2, resolution};
      auto idx = _block_index.find(key);
      if (idx != _block_index.end()) {
        idx->second.emplace(std::move(mm));
      } else {
        _block_index.emplace(key, phmap::btree_set<MatrixBlockMetadata>{std::move(mm)});
      }
    }
  } catch (...) {
    early_return = true;
    throw;
  }
}

}  // namespace hictk::hic::internal
