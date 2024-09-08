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
#include "hictk/filestream.hpp"
#include "hictk/hic.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/reference.hpp"
#include "hictk/transformers/coarsen.hpp"
#include "hictk/version.hpp"

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
  const auto match = _tank.find({chrom1, chrom2});
  if (match != _tank.end()) {
    return match->second;
  }
  throw std::out_of_range(
      fmt::format(FMT_STRING("MatrixBodyMetadataTank does not contain metadata for {}:{}"),
                  chrom1.name(), chrom2.name()));
}

inline HiCSectionOffsets MatrixBodyMetadataTank::offset(const Chromosome &chrom1,
                                                        const Chromosome &chrom2) const {
  const auto match = _offsets.find(Key{chrom1, chrom2});
  if (match != _offsets.end()) {
    return match->second;
  }
  throw std::out_of_range(
      fmt::format(FMT_STRING("MatrixBodyMetadataTank does not contain file offsets for {}:{}"),
                  chrom1.name(), chrom2.name()));
}

inline void MatrixBodyMetadataTank::insert(const Chromosome &chrom1, const Chromosome &chrom2,
                                           MatrixMetadata matrix_metadata,
                                           MatrixResolutionMetadata matrix_resolution_metadata) {
  try {
    auto match = _tank.find(Key{chrom1, chrom2});
    if (match != _tank.end()) {
      match->second.matrixMetadata = std::move(matrix_metadata);
      match->second.resolutionMetadata.emplace(std::move(matrix_resolution_metadata));
    } else {
      _tank.emplace(Key{chrom1, chrom2},
                    MatrixBodyMetadata{std::move(matrix_metadata),
                                       phmap::btree_set<MatrixResolutionMetadata>{
                                           std::move(matrix_resolution_metadata)}});
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while inserting metadata for {}:{} into a "
                               "MatrixBodyMetadataTank object: {}"),
                    chrom1.name(), chrom2.name(), e.what()));
  }
}

inline void MatrixBodyMetadataTank::update_offsets(const Chromosome &chrom1,
                                                   const Chromosome &chrom2,
                                                   std::streamoff position, std::size_t size) {
  try {
    auto [it, inserted] =
        _offsets.try_emplace(Key{chrom1, chrom2}, HiCSectionOffsets{position, size});
    if (!inserted) {
      it->second = HiCSectionOffsets{position, size};
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while updating MatrixBodyMetadata file offsets for {}:{}: {}"),
        chrom1.name(), chrom2.name(), e.what()));
  }
}

inline void MatrixBodyMetadataTank::remove(const Chromosome &chrom1, const Chromosome &chrom2) {
  try {
    const Key k{chrom1, chrom2};
    _tank.erase(k);
    _offsets.erase(k);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while removing the MatrixBodyMetadata entry for {}:{}: {}"),
        chrom1.name(), chrom2.name(), e.what()));
  }
}

inline auto MatrixBodyMetadataTank::operator()() const noexcept
    -> const phmap::flat_hash_map<Key, MatrixBodyMetadata> & {
  return _tank;
}

inline HiCFileWriter::HiCFileWriter(std::string_view path_, std::size_t n_threads)
    : _fs(std::string{path_}, std::ios::in | std::ios::out),
      _header(read_header(_fs)),
      _bin_tables(init_bin_tables(chromosomes(), resolutions())),
      _tpool(init_tpool(n_threads)) {
  read_offsets();
  read_norm_expected_values();
  read_norm_vectors();
}

inline HiCFileWriter::HiCFileWriter(std::string_view path_, Reference chromosomes_,
                                    std::vector<std::uint32_t> resolutions_,
                                    std::string_view assembly_, std::size_t n_threads,
                                    std::size_t chunk_size, const std::filesystem::path &tmpdir,
                                    std::uint32_t compression_lvl, bool skip_all_vs_all_matrix,
                                    std::size_t buffer_size)
    : _fs(filestream::FileStream::create(std::string{path_})),
      _tmpdir(tmpdir),
      _header(init_header(path_, std::move(chromosomes_), std::move(resolutions_), assembly_,
                          skip_all_vs_all_matrix)),
      _bin_tables(init_bin_tables(chromosomes(), resolutions())),
      _block_mappers(init_interaction_block_mappers(_tmpdir, _bin_tables, chunk_size, 3)),
      _compression_lvl(compression_lvl),
      _compressor(libdeflate_alloc_compressor(static_cast<std::int32_t>(compression_lvl))),
      _compression_buffer(buffer_size, '\0'),
      _tpool(init_tpool(n_threads)),
      _skip_all_vs_all_matrix(skip_all_vs_all_matrix) {
  if (!std::filesystem::exists(_tmpdir)) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("temporary directory {} does not exist"), _tmpdir));
  }
}

inline std::string_view HiCFileWriter::path() const noexcept { return _header.url; }

inline const Reference &HiCFileWriter::chromosomes() const noexcept { return _header.chromosomes; }

inline const BinTable &HiCFileWriter::bins(std::uint32_t resolution) const {
  return *_bin_tables.at(resolution);
}

inline const std::vector<std::uint32_t> &HiCFileWriter::resolutions() const noexcept {
  return _header.resolutions;
}

inline auto HiCFileWriter::stats(std::uint32_t resolution) const noexcept -> Stats {
  auto match = _stats.find(resolution);
  if (match != _stats.end()) {
    return match->second;
  }
  return {};
}

inline void HiCFileWriter::serialize() {
  try {
    write_header();
    write_pixels(_skip_all_vs_all_matrix);
    finalize(true);
    for (auto &[_, mapper] : _block_mappers) {
      mapper.clear();
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing file \"{}\": {}"), path(), e.what()));
  }
}

inline void HiCFileWriter::write_header() {
  assert(_fs.tellp() == 0);

  assert(_header.version == 9);
  assert(!chromosomes().empty());

  try {
    const auto offset1 = _fs.tellp();
    SPDLOG_INFO(FMT_STRING("writing header at offset {}"), offset1);
    _fs.write(_header.serialize(_bbuffer));
    const auto offset2 = _fs.tellp();

    _header_section = {offset1, offset2 - offset1};
    _data_block_section = {offset2, 0};
    _body_metadata_section = {offset2, 0};
    _footer_section = {offset2, 0};
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing the .hic header for file \"{}\" to disk: {}"),
        path(), e.what()));
  }
}

inline void HiCFileWriter::write_footer_size() {
  SPDLOG_DEBUG(FMT_STRING("updating footer size to {}"), _footer_section.size());
  // This is not documented for v9, but nBytesV5 is not included in the footer size
  const auto nBytesV5 = static_cast<std::int64_t>(_footer_section.size()) -
                        static_cast<std::int64_t>(sizeof(std::int64_t));

  try {
    _fs.seekp(_footer_section.start());
    _fs.write(nBytesV5);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing the footer size for file \"{}\" to disk: {}"),
        path(), e.what()));
  }
}

inline void HiCFileWriter::write_footer_offset() {
  SPDLOG_DEBUG(FMT_STRING("updating footer offset to {}"), _footer_section.start());
  const auto offset = sizeof("HIC") + sizeof(_header.version);

  try {
    _fs.seekp(offset);
    _fs.write(conditional_static_cast<std::int64_t>(_footer_section.start()));
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing the footer offset for file \"{}\" to disk: {}"),
        path(), e.what()));
  }
}

inline void HiCFileWriter::write_norm_vector_index() {
  const auto offset =
      static_cast<std::int64_t>(sizeof("HIC") + sizeof(_header.version) +
                                sizeof(_header.footerPosition) + _header.genomeID.size() + 1);
  const auto normVectorIndexPosition =
      conditional_static_cast<std::int64_t>(_norm_vector_index_section.start());
  const auto normVectorIndexLength = static_cast<std::int64_t>(_norm_vector_index_section.size());

  SPDLOG_DEBUG(FMT_STRING("writing normVectorIndex {}:{} at offset {}..."), normVectorIndexPosition,
               normVectorIndexLength, offset);

  try {
    _fs.seekp(offset);
    _fs.write(normVectorIndexPosition);
    _fs.write(normVectorIndexLength);
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while writing the normVectorIndex position and "
                               "length for file \"{}\" to disk: {}"),
                    path(), e.what()));
  }
}

template <typename PixelIt, typename>
inline void HiCFileWriter::add_pixels(std::uint32_t resolution, PixelIt first_pixel,
                                      PixelIt last_pixel) {
  try {
    _block_mappers.at(resolution).append_pixels(first_pixel, last_pixel, _tpool);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while adding pixels for resolution {} to file \"{}\": {}"),
        resolution, path(), e.what()));
  }
}

inline void HiCFileWriter::write_pixels(bool skip_all_vs_all_matrix) {
  SPDLOG_INFO(FMT_STRING("begin writing interaction blocks to file \"{}\"..."), path());
  const auto &chrom_idx = _block_mappers.at(resolutions().front()).chromosome_index();
  std::vector<std::pair<Chromosome, Chromosome>> chroms{chrom_idx.size()};
  std::transform(chrom_idx.begin(), chrom_idx.end(), chroms.begin(),
                 [](const auto &kv) { return kv.first; });
  std::sort(chroms.begin(), chroms.end());

  for (const auto &[chrom1, chrom2] : chroms) {
    if (chrom1.is_all() || chrom2.is_all()) {
      continue;
    }
    write_pixels(chrom1, chrom2);
  }

  if (!skip_all_vs_all_matrix) {
    write_all_matrix();
  }
}

inline void HiCFileWriter::write_all_matrix(std::uint32_t target_num_bins) {
  try {
    std::uint64_t genome_size = 0;
    for (const auto &chrom : chromosomes()) {
      if (chrom.is_all()) {
        continue;
      }
      genome_size += chrom.size();
    }

    const auto base_resolution = resolutions().back();
    auto target_resolution =
        static_cast<std::uint32_t>((genome_size + target_num_bins - 1) / target_num_bins);
    const auto factor = std::max(std::uint32_t(1), target_resolution / base_resolution);
    target_resolution = factor * base_resolution;
    const auto target_resolution_scaled =
        std::max(std::uint32_t{1}, target_resolution / DEFAULT_CHROM_ALL_SCALE_FACTOR);

    SPDLOG_INFO(FMT_STRING("writing pixels for {}:{} matrix..."), chromosomes().at(0).name(),
                chromosomes().at(0).name());

    std::uint32_t genome_size_scaled = 0;
    for (const auto &chrom : chromosomes()) {
      if (chrom.is_all()) {
        continue;
      }
      const auto num_bins = (chrom.size() + target_resolution - 1) / target_resolution;
      genome_size_scaled += static_cast<std::uint32_t>(num_bins) * target_resolution_scaled;
    }

    genome_size_scaled = std::max(std::uint32_t{1}, genome_size_scaled);

    const auto bin_table_ALL = std::make_shared<const BinTable>(
        Reference{Chromosome{0, "__ALL__", genome_size_scaled}}, target_resolution_scaled);
    const auto chrom = bin_table_ALL->chromosomes().at(0);

    const auto num_bins =
        HiCInteractionToBlockMapper::compute_num_bins(chrom, chrom, target_resolution_scaled);
    const auto num_columns = HiCInteractionToBlockMapper::compute_block_column_count(
        chrom, chrom, target_resolution_scaled, HiCInteractionToBlockMapper::DEFAULT_INTER_CUTOFF);
    const auto num_rows = num_bins / num_columns + 1;

    HiCInteractionToBlockMapper::BlockMapperIntra mapper{num_rows, num_columns};

    const File f(std::string{path()}, base_resolution);
    auto sel = f.fetch();
    phmap::btree_map<std::uint64_t, MatrixInteractionBlock<float>> blocks{};

    std::for_each(sel.begin<float>(), sel.end<float>(), [&](const ThinPixel<float> &p) {
      const Pixel<float> pixel(*_bin_tables.at(base_resolution), p);
      // The result of this coarsening is not correct, as the last bin in a chromosome will
      // have the same ID as the first bin in the next chromosome, but this is what JuiceBox
      // expects.
      // We subtract the chromosome ID as JuiceBox's chromosome grid expects pixels boundaries to be
      // multiples of the bin size. This turns out to be correct as long as chromosome sizes are not
      // multiples of the bin size (which should happen extremely rarely), in which case the result
      // is off by one.
      Pixel<float> coarsened_pixel(
          *bin_table_ALL, (p.bin1_id - (pixel.coords.bin1.chrom().id() - 1)) / factor,
          (p.bin2_id - (pixel.coords.bin2.chrom().id() - 1)) / factor, p.count);

      const auto bid =
          mapper(coarsened_pixel.coords.bin1.rel_id(), coarsened_pixel.coords.bin2.rel_id());
      auto [it, inserted] = blocks.try_emplace(bid, MatrixInteractionBlock<float>{});
      it->second.emplace_back(std::move(coarsened_pixel));
    });

    const auto offset = _data_block_section.end();
    _fs.seekp(offset);

    for (auto &[bid, blk] : blocks) {
      blk.finalize();
      write_interaction_block(bid, chrom, chrom, target_resolution_scaled, blk);
    }
    _data_block_section.size() += _fs.tellp() - static_cast<std::size_t>(offset);

    add_body_metadata(target_resolution_scaled, chrom, chrom);
    write_body_metadata();
    add_footer(chrom, chrom);
    write_footers();

    finalize();
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing the All:All matrix to file \"{}\": {}"), path(),
        e.what()));
  }
}

inline auto HiCFileWriter::write_pixels(const Chromosome &chrom1,
                                        const Chromosome &chrom2) -> HiCSectionOffsets {
  try {
    write_pixels(chrom1, chrom2, resolutions().front());
    add_body_metadata(resolutions().front(), chrom1, chrom2);
    write_body_metadata();
    add_footer(chrom1, chrom2);
    write_footers();

    finalize();
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while writing the {}:{} matrix at {} resolution to file \"{}\": {}"),
        chrom1.name(), chrom2.name(), resolutions().front(), path(), e.what()));
  }

  for (std::size_t i = 1; i < resolutions().size(); ++i) {
    auto base_resolution = resolutions().front();
    const auto res = resolutions()[i];

    auto &mapper = _block_mappers.at(res);
    if (mapper.empty(chrom1, chrom2)) {
      try {
        for (std::size_t j = 0; j < i; ++j) {
          if (res % resolutions()[j] == 0) {
            base_resolution = resolutions()[j];
          }
        }
        const File f(std::string{path()}, base_resolution);
        const auto sel = f.fetch(chrom1.name(), chrom2.name());
        if (!sel.empty()) {
          SPDLOG_INFO(
              FMT_STRING("[{} bp] no pixels provided for {}:{} matrix: generating pixels by "
                         "coarsening resolution {}..."),
              res, chrom1.name(), chrom2.name(), base_resolution);
          const auto factor = res / base_resolution;
          const transformers::CoarsenPixels coarsener(
              sel.begin<float>(), sel.end<float>(),
              std::make_shared<const BinTable>(bins(base_resolution)), factor);

          mapper.append_pixels(coarsener.begin(), coarsener.end(), _tpool);
        }
      } catch (const std::exception &e) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("an error occurred while coarsening interactions for {}:{} from "
                                   "resolution {} to resolution {}: {}"),
                        chrom1.name(), chrom2.name(), base_resolution, res, e.what()));
      }
    }

    if (mapper.empty(chrom1, chrom2)) {
      SPDLOG_WARN(FMT_STRING("[{} bp] no pixels found for {}:{} matrix: SKIPPING!"), res,
                  chrom1.name(), chrom2.name());
      continue;
    }

    try {
      mapper.finalize();
      write_pixels(chrom1, chrom2, res);
      for (std::size_t j = 0; j <= i; ++j) {
        add_body_metadata(resolutions()[j], chrom1, chrom2);
      }
      write_body_metadata();
      add_footer(chrom1, chrom2);
      write_footers();
      finalize();
    } catch (const std::exception &e) {
      throw std::runtime_error(fmt::format(FMT_STRING("an error occurred while writing the {}:{} "
                                                      "matrix at {} resolution to file \"{}\": {}"),
                                           chrom1.name(), chrom2.name(), res, path(), e.what()));
    }
  }
  return {_data_block_section.start(),
          _fs.tellp() - static_cast<std::size_t>(_data_block_section.start())};
}

inline void HiCFileWriter::write_body_metadata() {
  const auto pos = _data_block_section.end();
  _fs.seekp(pos);
  for (const auto &[chroms, metadata] : _matrix_metadata()) {
    const auto &chrom1 = chroms.chrom1;
    const auto &chrom2 = chroms.chrom2;
    [[maybe_unused]] const auto &num_resolutions = metadata.resolutionMetadata.size();

    try {
      const auto pos1 = _fs.tellp();
      SPDLOG_DEBUG(FMT_STRING("writing MatrixBodyMetadata for {}:{} ({} resolutions) at offset {}"),
                   chrom1.name(), chrom2.name(), num_resolutions, pos1);
      _fs.write(metadata.serialize(_bbuffer));
      const auto pos2 = _fs.tellp();
      SPDLOG_DEBUG(FMT_STRING("updating MatrixBodyMetadata offset and size for {}:{} ({} "
                              "resolutions) to {} and {}"),
                   chrom1.name(), chrom2.name(), num_resolutions, pos1, pos2 - pos1);
      _matrix_metadata.update_offsets(chrom1, chrom2, static_cast<std::streamoff>(pos1),
                                      pos2 - pos1);
    } catch (const std::exception &e) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("an error occurred while writing the MatrixBodyMetadata for {}:{} "
                                 "to file \"{}\": {}"),
                      chrom1.name(), chrom2.name(), path(), e.what()));
    }
  }

  const auto size = _fs.tellp() - static_cast<std::size_t>(pos);
  _body_metadata_section = {pos, size};
}

inline void HiCFileWriter::add_body_metadata(std::uint32_t resolution, const Chromosome &chrom1,
                                             const Chromosome &chrom2, const std::string &unit) {
  SPDLOG_DEBUG(FMT_STRING("adding MatrixBodyMetadata for {}:{} at {} {}"), chrom1.name(),
               chrom2.name(), resolution, unit);
  const auto sum_counts =
      chrom1.name() == "__ALL__" ? 1.0F : _block_mappers.at(resolution).pixel_sum(chrom1, chrom2);
  if (sum_counts == 0) {
    return;
  }

  try {
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

    const auto &blks = _block_index.at(BlockIndexKey{chrom1, chrom2, resolution});
    mrm.set_block_metadata(blks.begin(), blks.end());

    mm.chr1Idx = static_cast<std::int32_t>(chrom1.id());
    mm.chr2Idx = static_cast<std::int32_t>(chrom2.id());
    mm.nResolutions = static_cast<std::int32_t>(metadata.resolutionMetadata.size() + 1);

    _matrix_metadata.insert(chrom1, chrom2, mm, mrm);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while adding or updating the MatrixBodyMetadata for {}:{}: {}"),
        chrom1.name(), chrom2.name(), e.what()));
  }
}

inline void HiCFileWriter::write_footers() {
  const auto offset1 = _body_metadata_section.end();

  try {
    _fs.seekp(offset1);
    SPDLOG_DEBUG(FMT_STRING("initializing footer section at offset {}"), offset1);
    const std::int64_t nBytesV5 = -1;
    const auto nEntries = static_cast<std::int32_t>(_footers.size());
    _fs.write(nBytesV5);
    _fs.write(nEntries);

    for (auto &[chroms, footer] : _footers) {
      try {
        const auto offset = _matrix_metadata.offset(chroms.first, chroms.second);
        footer.position = conditional_static_cast<std::int64_t>(offset.start());
        footer.size = static_cast<std::int32_t>(offset.size());
        SPDLOG_DEBUG(FMT_STRING("writing FooterMasterIndex for {}:{} at offset {}"),
                     chroms.first.name(), chroms.second.name(), _fs.tellp());
        _fs.write(footer.serialize(_bbuffer));
      } catch (const std::exception &e) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("an error occurred while writing the footer for {}:{}: {}"),
                        chroms.first.name(), chroms.second.name(), e.what()));
      }
    }

    write_empty_expected_values();
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing the footer section to file \"{}\": {}"), path(),
        e.what()));
  }

  _footer_section = {offset1, _fs.tellp() - static_cast<std::size_t>(offset1)};
}

inline void HiCFileWriter::add_footer(const Chromosome &chrom1, const Chromosome &chrom2) {
  if (!_matrix_metadata.contains(chrom1, chrom2)) {
    return;
  }

  try {
    FooterMasterIndex footer{};
    footer.key = fmt::format(FMT_STRING("{}_{}"), chrom1.id(), chrom2.id());
    footer.position = -1;
    footer.size = -1;

    auto [it, inserted] = _footers.emplace(std::make_pair(chrom1, chrom2), footer);
    if (!inserted) {
      it->second = std::move(footer);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while adding the footer for {}:{}: {}"),
                    chrom1.name(), chrom2.name(), e.what()));
  }
}

inline void HiCFileWriter::write_norm_vectors_and_norm_expected_values() {
  // we are writing the norm vectors twice because the function computing the norm expected values
  // expects the normalization vectors to be available in the file that is being written
  write_norm_vectors();
  compute_and_write_normalized_expected_values();
  write_norm_vectors();
}

inline void HiCFileWriter::write_empty_expected_values() {
  ExpectedValues ev{};

  try {
    const auto offset = _fs.tellp();
    SPDLOG_DEBUG(FMT_STRING("writing empty expected values section at offset {}..."), offset);
    _fs.write(ev.serialize(_bbuffer));

    _expected_values_section = {offset, _fs.tellp() - offset};
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while writing an empty expected values section to file \"{}\": {}"),
        path(), e.what()));
  }
}

inline void HiCFileWriter::write_empty_normalized_expected_values() {
  const auto offset = _expected_values_section.end();
  SPDLOG_DEBUG(FMT_STRING("writing empty expected values (normalized) section at offset {}..."),
               offset);
  try {
    _fs.seekp(offset);
    HICTK_DISABLE_WARNING_PUSH
    HICTK_DISABLE_WARNING_USELESS_CAST
    _fs.write(std::int32_t(0));
    HICTK_DISABLE_WARNING_POP
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while writing an empty normalized expected "
                               "values section to file \"{}\": {}"),
                    path(), e.what()));
  }
  _expected_values_norm_section = {offset, _fs.tellp() - static_cast<std::size_t>(offset)};
}

inline ExpectedValuesBlock HiCFileWriter::compute_expected_values(std::uint32_t resolution) {
  SPDLOG_DEBUG(FMT_STRING("computing expected values at resolution {}..."), resolution);

  try {
    const File f(std::string{path()}, resolution);
    const auto sel = f.fetch();

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

    return {"BP", resolution, aggr.weights(), chrom_ids, scaling_factors};
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while computing the expected values for file "
                               "\"{}\" at {} resolution: {}"),
                    path(), resolution, e.what()));
  }
}

inline NormalizedExpectedValuesBlock HiCFileWriter::compute_normalized_expected_values(
    std::uint32_t resolution, const balancing::Method &norm) {
  assert(norm != balancing::Method::NONE());
  SPDLOG_INFO(FMT_STRING("computing normalized expected values ({}) at resolution {}..."), norm,
              resolution);

  try {
    const File f(std::string{path()}, resolution);
    const auto sel = f.fetch(norm);

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

    return {norm.to_string(), "BP", resolution, aggr.weights(), chrom_ids, scaling_factors};
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while computing the normalized expected values for file "
                   "\"{}\" at {} resolution: {}"),
        path(), resolution, e.what()));
  }
}

inline void HiCFileWriter::compute_and_write_expected_values() {
  assert(_tpool.get_thread_count() != 0);
  ExpectedValues ev{};

  std::vector<std::future<ExpectedValuesBlock>> results{};
  for (const auto &resolution : resolutions()) {
    results.emplace_back(
        _tpool.submit_task([&, res = resolution]() { return compute_expected_values(res); }));
  }

  for (auto &res : results) {
    ev.emplace(res.get());
  }

  try {
    const auto offset = _footer_section.end() -
                        conditional_static_cast<std::streamoff>(sizeof(ev.nExpectedValueVectors()));
    SPDLOG_INFO(FMT_STRING("writing {} expected value vectors at offset {}..."),
                ev.nExpectedValueVectors(), offset);
    _fs.seekp(offset);
    _fs.write(ev.serialize(_bbuffer));

    _expected_values_section = {offset, _fs.tellp() - static_cast<std::size_t>(offset)};
    _footer_section.size() += _expected_values_section.size() - sizeof(ev.nExpectedValueVectors());
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing expected values to file \"{}\": {}"), path(),
        e.what()));
  }
}

inline void HiCFileWriter::compute_and_write_normalized_expected_values() {
  assert(_tpool.get_thread_count() != 0);
  NormalizedExpectedValues ev{};

  phmap::btree_map<NormalizedExpectedValuesBlock, std::future<NormalizedExpectedValuesBlock>>
      results{};
  for (const auto &[blk, _] : _normalization_vectors) {
    const NormalizedExpectedValuesBlock key{
        blk.type, blk.unit, static_cast<std::uint32_t>(blk.binSize), {}, {}, {}};
    const auto nev_available =
        _normalized_expected_values.find(key) != _normalized_expected_values.end();
    const auto nev_already_submitted_for_computation = results.find(key) != results.end();
    if (!nev_available && !nev_already_submitted_for_computation) {
      results.emplace(key, _tpool.submit_task([&, res = static_cast<std::uint32_t>(blk.binSize),
                                               type = blk.type]() {
        const balancing::Method norm{type};
        return compute_normalized_expected_values(res, norm);
      }));
    }
  }

  for (auto &[_, res] : results) {
    _normalized_expected_values.emplace(res.get());
  }

  for (const auto &nev : _normalized_expected_values) {
    ev.emplace(nev);
  }

  try {
    const auto offset = _footer_section.end();
    SPDLOG_INFO(FMT_STRING("writing {} normalized expected value vectors at offset {}..."),
                ev.nNormExpectedValueVectors(), offset);
    _fs.seekp(offset);
    _fs.write(ev.serialize(_bbuffer));

    _expected_values_norm_section = {offset, _fs.tellp() - static_cast<std::size_t>(offset)};
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing normalized expected values to file \"{}\": {}"),
        path(), e.what()));
  }
}

inline void HiCFileWriter::add_norm_vector(const NormalizationVectorIndexBlock &blk,
                                           const std::vector<float> &weights,
                                           bool force_overwrite) {
  const auto &chrom = chromosomes().at(static_cast<std::uint32_t>(blk.chrIdx));
  SPDLOG_INFO(FMT_STRING("[{}] adding {} normalization vector for {} ({}): {} values"), blk.binSize,
              blk.type, chrom.name(), blk.unit, weights.size());

  try {
    const auto bin_size = static_cast<std::uint32_t>(blk.binSize);
    const auto expected_shape = (chrom.size() + bin_size - 1) / bin_size;

    if (weights.size() != expected_shape) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("weight shape mismatch: expected {} values, found {}"),
                      expected_shape, weights.size()));
    }

    auto [it, inserted] = _normalization_vectors.emplace(blk, weights);
    if (!inserted) {
      if (force_overwrite) {
        it->second = weights;
        const NormalizedExpectedValuesBlock key{
            blk.type, blk.unit, static_cast<std::uint32_t>(blk.binSize), {}, {}, {}};
        _normalized_expected_values.erase(key);
      } else {
        throw std::runtime_error("file already contains normalization vector");
      }
    }

  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while adding {} normalization vector for {} at {} resolution: {}"),
        blk.type, chrom.name(), blk.binSize, e.what()));
  }
}

inline void HiCFileWriter::add_norm_vector(std::string_view type, const Chromosome &chrom,
                                           std::string_view unit, std::uint32_t bin_size,
                                           const balancing::Weights &weights, bool force_overwrite,
                                           std::size_t position, std::size_t n_bytes) {
  add_norm_vector(NormalizationVectorIndexBlock{std::string{type}, chrom.id(), std::string{unit},
                                                bin_size, position, n_bytes},
                  weights, force_overwrite);
}

inline void HiCFileWriter::add_norm_vector(const NormalizationVectorIndexBlock &blk,
                                           const balancing::Weights &weights,
                                           bool force_overwrite) {
  std::vector<float> weights_f(weights.size());
  const auto weights_d = weights(balancing::Weights::Type::DIVISIVE);
  std::transform(weights_d.begin(), weights_d.end(), weights_f.begin(),
                 [](const auto n) { return static_cast<float>(n); });
  add_norm_vector(blk, weights_f, force_overwrite);
}

inline void HiCFileWriter::add_norm_vector(std::string_view type, std::string_view unit,
                                           std::uint32_t bin_size,
                                           const balancing::Weights &weights,
                                           bool force_overwrite) {
  try {
    const auto expected_shape = bins(bin_size).size();
    if (weights.size() != expected_shape) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("weight shape mismatch: expected {} values, found {}"),
                      expected_shape, weights.size()));
    }

    std::ptrdiff_t i0 = 0;
    std::ptrdiff_t i1 = 0;
    const auto weights_ = weights(balancing::Weights::Type::DIVISIVE);
    for (const auto &chrom : chromosomes()) {
      if (chrom.is_all()) {
        continue;
      }
      i1 += static_cast<std::ptrdiff_t>((chrom.size() + bin_size - 1) / bin_size);
      std::vector<double> chrom_weights(weights_.begin() + i0, weights_.begin() + i1);
      add_norm_vector(type, chrom, unit, bin_size,
                      {std::move(chrom_weights), balancing::Weights::Type::DIVISIVE},
                      force_overwrite);
      i0 = i1;
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("an error occurred while adding {} genome-wide "
                                                    "normalization vector at {} resolution: {}"),
                                         type, bin_size, e.what()));
  }
}

inline void HiCFileWriter::finalize(bool compute_expected_values) {
  try {
    if (compute_expected_values) {
      compute_and_write_expected_values();
      write_empty_normalized_expected_values();
      write_norm_vectors();
      compute_and_write_normalized_expected_values();
    } else {
      write_empty_expected_values();
      write_empty_normalized_expected_values();
    }

    write_footer_offset();
    write_footer_size();
    write_norm_vectors();
    _fs.flush();
    _fs.seekp(0, std::ios::end);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while finalizing file \"{}\": {}"), path(), e.what()));
  }
}

inline void HiCFileWriter::write_norm_vectors() {
  try {
    const auto offset1 =
        std::max(_expected_values_norm_section.end(), _norm_vector_index_section.start());
    _fs.seekp(offset1);

    if (_normalization_vectors.empty()) {
      SPDLOG_DEBUG(FMT_STRING("writing empty normalization vector section at offset {}..."),
                   offset1);
    } else {
      SPDLOG_INFO(FMT_STRING("writing {} normalization vectors at offset {}..."),
                  _normalization_vectors.size(), offset1);
    }

    const auto nNormVectors = static_cast<std::int32_t>(_normalization_vectors.size());
    _fs.write(nNormVectors);

    phmap::btree_map<NormalizationVectorIndexBlock, HiCSectionOffsets> index_offsets{};
    for (const auto &[blk, _] : _normalization_vectors) {
      try {
        const auto offset2 = _fs.tellp();
        _fs.write(blk.serialize(_bbuffer));
        index_offsets.emplace(blk, HiCSectionOffsets{offset2, _fs.tellp() - offset2});
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(
            FMT_STRING(
                "an error occurred while writing the {} NormalizationVectorIndexBlock for {} "
                "at {} resolution to file \"{}\": {}"),
            blk.type, chromosomes().at(static_cast<std::uint32_t>(blk.chrIdx)).name(), blk.binSize,
            path(), e.what()));
      }
    }
    const auto offset2 = _fs.tellp();

    phmap::btree_map<NormalizationVectorIndexBlock, HiCSectionOffsets> vector_offsets{};
    for (const auto &[blk, weights] : _normalization_vectors) {
      try {
        const auto offset3 = _fs.tellp();
        const auto nValues = static_cast<std::int64_t>(weights.size());
        _fs.write(nValues);
        _fs.write(weights);
        vector_offsets.emplace(blk, HiCSectionOffsets{offset3, _fs.tellp() - offset3});
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("an error occurred while writing the {} normalization vector for {} "
                       "at {} resolution to file \"{}\": {}"),
            blk.type, chromosomes().at(static_cast<std::uint32_t>(blk.chrIdx)).name(), blk.binSize,
            path(), e.what()));
      }
    }

    const auto offset4 = _fs.tellp();

    for (const auto &[blk, idx_offsets] : index_offsets) {
      try {
        const auto &vect_offsets = vector_offsets.at(blk);
        auto new_blk = blk;
        new_blk.position = vect_offsets.start();
        new_blk.nBytes = static_cast<std::int64_t>(vect_offsets.size());
        _fs.seekp(idx_offsets.start());
        _fs.write(new_blk.serialize(_bbuffer));
      } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("an error occurred while updating file offsets in the {} "
                       "NormalizationVectorIndexBlock for {} at {} resolution to file \"{}\": {}"),
            blk.type, chromosomes().at(static_cast<std::uint32_t>(blk.chrIdx)).name(), blk.binSize,
            path(), e.what()));
      }
    }

    _norm_vector_index_section = {offset1, offset2 - static_cast<std::size_t>(offset1)};
    _norm_vectors_section = {offset2, offset4 - static_cast<std::size_t>(offset2)};

    write_norm_vector_index();
    _fs.seekp(0, std::ios::end);
    _fs.flush();
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing normalization vectors to file \"{}\": {}"),
        path(), e.what()));
  }
}

inline HiCHeader HiCFileWriter::read_header(filestream::FileStream &fs) {
  return HiCHeader::deserialize(fs);
}

inline HiCHeader HiCFileWriter::init_header(std::string_view path, Reference chromosomes,
                                            std::vector<std::uint32_t> resolutions,
                                            std::string_view assembly,
                                            bool skip_all_vs_all_matrix) {
  if (skip_all_vs_all_matrix) {
    chromosomes = chromosomes.remove_ALL();
  } else {
    chromosomes = chromosomes.add_ALL(DEFAULT_CHROM_ALL_SCALE_FACTOR);
  }
  return {
      std::string{path},      // url
      9,                      // version
      -1,                     // footerPosition
      std::string{assembly},  // genomeId
      -1,                     // normVectorIndexPosition
      0,                      // normVectorIndexLength
      std::move(chromosomes),
      std::move(resolutions),                                   // resolutions
      {{"software", std::string{config::version::str_long()}}}  // attributes
  };
}

inline auto HiCFileWriter::init_bin_tables(
    const Reference &chromosomes, const std::vector<std::uint32_t> &resolutions) -> BinTables {
  BinTables bin_tables(resolutions.size());
  for (const auto &res : resolutions) {
    bin_tables.emplace(res, std::make_shared<const BinTable>(chromosomes, res));
  }

  return bin_tables;
}

inline auto HiCFileWriter::init_interaction_block_mappers(const std::filesystem::path &root_folder,
                                                          const BinTables &bin_tables,
                                                          std::size_t chunk_size,
                                                          int compression_lvl) -> BlockMappers {
  BlockMappers mappers(bin_tables.size());
  for (const auto &[res, bin_table] : bin_tables) {
    const auto path = fmt::format(FMT_STRING("{}/{}.bin"), root_folder.string(), res);
    mappers.emplace(res, HiCInteractionToBlockMapper{path, bin_table, chunk_size, compression_lvl});
  }

  return mappers;
}

inline BS::thread_pool HiCFileWriter::init_tpool(std::size_t n_threads) {
  return BS::thread_pool{
      conditional_static_cast<BS::concurrency_t>(n_threads < 2 ? std::size_t(1) : n_threads)};
}

inline auto HiCFileWriter::write_pixels(const Chromosome &chrom1, const Chromosome &chrom2,
                                        std::uint32_t resolution) -> HiCSectionOffsets {
  try {
    const auto offset = _data_block_section.end();
    _fs.seekp(offset);

    SPDLOG_INFO(FMT_STRING("[{} bp] writing pixels for {}:{} matrix at offset {}..."), resolution,
                chrom1.name(), chrom2.name(), offset);

    const auto stats = write_interaction_blocks(chrom1, chrom2, resolution);

    SPDLOG_INFO(FMT_STRING("[{} bp] written {} pixels for {}:{} matrix"), resolution, stats.nnz,
                chrom1.name(), chrom2.name());

    auto [it, inserted] = _stats.try_emplace(resolution, stats);
    if (!inserted) {
      it->second.sum += stats.sum;
      it->second.nnz += stats.nnz;
    }

    _data_block_section.size() += _fs.tellp() - static_cast<std::size_t>(offset);
    return {offset, _fs.tellp() - static_cast<std::size_t>(offset)};
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while writing pixels for {}:{} to file \"{}\": {}"),
        chrom1.name(), chrom2.name(), path(), e.what()));
  }
}

inline auto HiCFileWriter::write_interaction_blocks(const Chromosome &chrom1,
                                                    const Chromosome &chrom2,
                                                    std::uint32_t resolution) -> Stats {
  auto &mapper = _block_mappers.at(resolution);
  mapper.finalize();

  const auto block_ids = mapper.chromosome_index().find(std::make_pair(chrom1, chrom2));
  if (block_ids == mapper.chromosome_index().end()) {
    SPDLOG_DEBUG(FMT_STRING("no pixels to write for {}:{} matrix at {} resolution"), chrom1.name(),
                 chrom2.name(), resolution);
    return {};
  }

  if (_tpool.get_thread_count() < 3 || block_ids->second.size() == 1) {
    try {
      Stats stats{};
      for (const auto &bid : block_ids->second) {
        auto blk = mapper.merge_blocks(bid);
        stats.sum += blk.sum();
        stats.nnz += blk.size();
        write_interaction_block(bid.bid, chrom1, chrom2, resolution, std::move(blk));
      }

      return stats;
    } catch (const std::exception &e) {
      throw std::runtime_error(
          "an error occurred while writing interaction blocks using a single thread: " +
          std::string{e.what()});
    }
  }

  try {
    std::mutex block_id_queue_mtx{};
    std::queue<std::uint64_t> block_id_queue{};
    moodycamel::BlockingConcurrentQueue<HiCInteractionToBlockMapper::BlockID> block_queue(
        block_ids->second.size());

    std::mutex serialized_block_tank_mtx{};
    phmap::flat_hash_map<std::uint64_t, std::string> serialized_block_tank{
        block_ids->second.size()};
    const auto stop_token = std::numeric_limits<std::uint64_t>::max();
    std::atomic<bool> early_return = false;

    std::mutex mapper_mtx{};

    std::vector<std::future<Stats>> worker_threads{};
    for (BS::concurrency_t i = 2; i < _tpool.get_thread_count(); ++i) {
      worker_threads.emplace_back(_tpool.submit_task([&]() {
        return merge_and_compress_blocks_thr(mapper, mapper_mtx, block_id_queue, block_id_queue_mtx,
                                             block_queue, serialized_block_tank,
                                             serialized_block_tank_mtx, early_return, stop_token);
      }));
    }

    auto writer = _tpool.submit_task([&]() {
      write_compressed_blocks_thr(chrom1, chrom2, resolution, block_id_queue, block_id_queue_mtx,
                                  serialized_block_tank, serialized_block_tank_mtx, early_return,
                                  stop_token);
    });

    auto producer = _tpool.submit_task([&]() {
      try {
        for (const auto &bid : block_ids->second) {
          if (early_return) {
            break;
          }

          block_queue.enqueue(bid);
          if (early_return) {
            break;
          }
        }
        for (std::size_t i = 0; i < worker_threads.size(); ++i) {
          block_queue.enqueue(HiCInteractionToBlockMapper::BlockID{0, 0, stop_token});
        }
      } catch (const std::exception &e) {
        early_return = true;
        throw std::runtime_error("an error occurred in the producer thread: " +
                                 std::string{e.what()});
      } catch (...) {
        early_return = true;
        throw;
      }
    });

    producer.get();

    Stats stats{};
    for (auto &worker : worker_threads) {
      const auto partial_stats = worker.get();
      stats.sum += partial_stats.sum;
      stats.nnz += partial_stats.nnz;
    }
    // signal no more blocks will be enqueued
    {
      std::scoped_lock lck(block_id_queue_mtx);
      block_id_queue.emplace(stop_token);
    }
    writer.get();

    return stats;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("an error occurred while interaction blocks using {} threads: {}"),
                    _tpool.get_thread_count(), e.what()));
  }
}

inline auto HiCFileWriter::write_interaction_block(
    std::uint64_t block_id, const Chromosome &chrom1, const Chromosome &chrom2,
    std::uint32_t resolution, const MatrixInteractionBlock<float> &blk) -> HiCSectionOffsets {
  const auto offset = _fs.tellp();

  std::ignore = blk.serialize(_bbuffer, *_compressor, _compression_buffer);
  SPDLOG_DEBUG(FMT_STRING("writing block #{} for {}:{}:{} at {}:{}"), block_id, chrom1.name(),
               chrom2.name(), resolution, offset, _compression_buffer.size());
  _fs.write(_compression_buffer);

  MatrixBlockMetadata mm{static_cast<std::int32_t>(block_id), static_cast<std::int64_t>(offset),
                         static_cast<std::int32_t>(_fs.tellp() - offset)};

  const BlockIndexKey key{chrom1, chrom2, resolution};
  auto idx = _block_index.find(key);
  if (idx != _block_index.end()) {
    idx->second.emplace(std::move(mm));
  } else {
    _block_index.emplace(key, phmap::btree_set<MatrixBlockMetadata>{std::move(mm)});
  }
  return {offset, _fs.tellp() - offset};
}

inline std::size_t HiCFileWriter::compute_num_bins(const Chromosome &chrom1,
                                                   const Chromosome &chrom2,
                                                   std::uint32_t resolution) {
  return HiCInteractionToBlockMapper::compute_num_bins(chrom1, chrom2, resolution);
}

inline void HiCFileWriter::add_norm_expected_values(const NormalizedExpectedValuesBlock &blk,
                                                    bool force_overwrite) {
  try {
    auto [it, inserted] = _normalized_expected_values.emplace(blk);
    if (!inserted) {
      if (force_overwrite) {
        *it = blk;
      } else {
        throw std::runtime_error("file already contains normalized expected values");
      }
    }

  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING(
            "an error occurred while adding {} normalized expected values at {} resolution: {}"),
        blk.type, blk.binSize, e.what()));
  }
}

inline void HiCFileWriter::read_norm_expected_values() {
  assert(_expected_values_norm_section.start() != 0);
  try {
    const auto offset = _expected_values_norm_section.start();
    _fs.seekg(offset);
    const auto nev = NormalizedExpectedValues::deserialize(_fs);

    for (const auto &ev : nev.normExpectedValues()) {
      add_norm_expected_values(ev);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(FMT_STRING("an error occurred while reading normalized "
                                                    "expected value vectors from file \"{}\": {}"),
                                         path(), e.what()));
  }
}

inline void HiCFileWriter::read_norm_vectors() {
  assert(_norm_vector_index_section.start() != 0);
  try {
    const auto offset = _norm_vector_index_section.start();
    _fs.seekg(offset);
    const auto nvi = NormalizationVectorIndex::deserialize(_fs);

    for (const auto &blk : nvi.normalizationVectorIndex()) {
      add_norm_vector(blk, read_norm_vector(blk), true);
    }
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while reading normalization vectors from file \"{}\": {}"),
        path(), e.what()));
  }
}

inline std::vector<float> HiCFileWriter::read_norm_vector(
    const NormalizationVectorIndexBlock &blk) {
  try {
    const auto offset = blk.position;
    _fs.seekg(offset);

    const auto &chrom = chromosomes().at(static_cast<std::uint32_t>(blk.chrIdx));
    const auto bin_size = static_cast<std::size_t>(blk.binSize);
    const auto nValuesExpected = (static_cast<std::size_t>(chrom.size()) + bin_size - 1) / bin_size;

    // https://github.com/aidenlab/hic-format/blob/master/HiCFormatV9.md#normalization-vector-arrays-1-per-normalization-vector
    const auto nValues = static_cast<std::size_t>(_fs.read<std::int64_t>());
    // We cannot use numValues directly because sometimes hic files have few trailing zeros for some
    // reason
    if (nValues < nValuesExpected) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("expected {} values, found {}"), nValuesExpected, nValues));
    }

    std::vector<float> buffer(nValues);
    _fs.read(buffer);
    buffer.resize(nValuesExpected);
    const auto bytes_read = _fs.tellg() - static_cast<std::size_t>(offset);
    if (bytes_read != static_cast<std::size_t>(blk.nBytes)) {
      throw std::runtime_error(
          fmt::format(FMT_STRING("expected to read {} bytes but read {}"), blk.nBytes, bytes_read));
    }
    return buffer;
  } catch (const std::exception &e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("{} normalization vector for {} at {} resolution is corrupted: {}"),
                    blk.type, _header.chromosomes.at(static_cast<std::uint32_t>(blk.chrIdx)).name(),
                    blk.binSize, e.what()));
  }
}

inline void HiCFileWriter::read_offsets() {
  try {
    _fs.seekg(0, std::ios::beg);
    const auto header_start = _fs.tellg();
    const auto header = HiCHeader::deserialize(_fs);
    const auto header_end = _fs.tellg();

    // read footer offsets
    _fs.seekg(header.footerPosition);
    const auto footer_start = _fs.tellg();
    const auto nBytesV5 = _fs.read<std::int64_t>();
    _fs.seekg(nBytesV5, std::ios::cur);
    const auto footer_end = _fs.tellg();

    // read norm expected values offsets
    const auto norm_expected_values_start = _fs.tellg();
    const auto nNormExpectedValueVectors = _fs.read<std::int32_t>();
    for (std::int32_t i = 0; i < nNormExpectedValueVectors; ++i) {
      std::ignore = NormalizationVectorIndexBlock::deserialize(_fs);
    }
    const auto norm_expected_values_end = _fs.tellg();

    // compute norm vector index offsets
    const auto norm_vector_index_start = header.normVectorIndexPosition;
    const auto norm_vector_index_end =
        header.normVectorIndexPosition + header.normVectorIndexLength;

    // set the offsets
    _header_section = {header_start, header_end - header_start};
    _footer_section = {footer_start, footer_end - footer_start};
    _expected_values_norm_section = {norm_expected_values_start,
                                     norm_expected_values_end - norm_expected_values_start};
    _norm_vector_index_section = {norm_vector_index_start,
                                  norm_vector_index_end - norm_vector_index_start};

    _fs.seekg(0, std::ios::end);
  } catch (const std::exception &e) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while reading section offsets from file \"{}\": {}"), path(),
        e.what()));
  }
}

inline std::size_t HiCFileWriter::compute_block_column_count(const Chromosome &chrom1,
                                                             const Chromosome &chrom2,
                                                             std::uint32_t resolution) {
  return HiCInteractionToBlockMapper::compute_block_column_count(
      chrom1, chrom2, resolution,
      chrom1 == chrom2 ? HiCInteractionToBlockMapper::DEFAULT_INTRA_CUTOFF
                       : HiCInteractionToBlockMapper::DEFAULT_INTER_CUTOFF);
}

inline auto HiCFileWriter::merge_and_compress_blocks_thr(
    HiCInteractionToBlockMapper &mapper, std::mutex &mapper_mtx,
    std::queue<std::uint64_t> &block_id_queue, std::mutex &block_id_queue_mtx,
    moodycamel::BlockingConcurrentQueue<HiCInteractionToBlockMapper::BlockID> &block_queue,
    phmap::flat_hash_map<std::uint64_t, std::string> &serialized_block_tank,
    std::mutex &serialized_block_tank_mtx, std::atomic<bool> &early_return,
    std::uint64_t stop_token) -> Stats {
  SPDLOG_DEBUG(FMT_STRING("merge_and_compress_blocks thread: start-up..."));
  try {
    HiCInteractionToBlockMapper::BlockID buffer{};
    BinaryBuffer bbuffer{};
    std::string compression_buffer(16'000'000, '\0');
    std::unique_ptr<libdeflate_compressor> libdeflate_compressor(
        libdeflate_alloc_compressor(static_cast<std::int32_t>(_compression_lvl)));
    std::unique_ptr<ZSTD_DCtx_s> zstd_dctx{ZSTD_createDCtx()};

    Stats stats{};
    while (!early_return) {
      // dequeue block
      if (!block_queue.wait_dequeue_timed(buffer, std::chrono::milliseconds(500))) {
        continue;
      }
      if (buffer.bid == stop_token) {
        SPDLOG_DEBUG(
            FMT_STRING("merge_and_compress_blocks thread: processed all blocks. Returning!"));
        return stats;
      }

      SPDLOG_DEBUG(
          FMT_STRING("merge_and_compress_blocks thread: merging partial blocks for block #{}"),
          buffer.bid);
      // read and merge partial blocks
      auto blk = mapper.merge_blocks(buffer, bbuffer, *zstd_dctx, compression_buffer, mapper_mtx);
      stats.nnz += blk.size();
      stats.sum += blk.sum();

      // compress and serialize block
      std::ignore = blk.serialize(bbuffer, *libdeflate_compressor, compression_buffer);

      // enqueue serialized block
      std::scoped_lock lck(serialized_block_tank_mtx, block_id_queue_mtx);
      SPDLOG_DEBUG(FMT_STRING("merge_and_compress_blocks thread: done processing block #{}"),
                   buffer.bid);
      serialized_block_tank.emplace(std::make_pair(buffer.bid, compression_buffer));
      block_id_queue.emplace(buffer.bid);
    }

    return stats;
  } catch (const std::exception &e) {
    early_return = true;
    throw std::runtime_error("an error occurred in merge_and_compress_blocks thread: " +
                             std::string{e.what()});

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
  SPDLOG_DEBUG(FMT_STRING("write_compressed_blocks thread: start-up..."));
  try {
    std::string buffer;

    while (!early_return) {
      const auto do_sleep = [&]() {
        std::scoped_lock lck(block_id_queue_mtx);
        return block_id_queue.empty();
      }();

      if (do_sleep) {
        SPDLOG_DEBUG(
            FMT_STRING("write_compressed_blocks thread: no blocks to consume. Sleeping..."));
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
        SPDLOG_DEBUG(FMT_STRING(
            "write_compressed_blocks thread: no more blocks to be processed. Returning!"));
        return;
      }

      SPDLOG_DEBUG(FMT_STRING("write_compressed_blocks thread: waiting for block #{}..."), bid);
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

      const auto offset = _fs.tellp();
      SPDLOG_DEBUG(FMT_STRING("writing block #{} for {}:{}:{} at {}:{}"), bid, chrom1.name(),
                   chrom2.name(), resolution, offset, buffer.size());
      _fs.write(buffer);

      MatrixBlockMetadata mm{static_cast<std::int32_t>(bid), static_cast<std::int64_t>(offset),
                             static_cast<std::int32_t>(buffer.size())};
      const BlockIndexKey key{chrom1, chrom2, resolution};
      auto idx = _block_index.find(key);
      if (idx != _block_index.end()) {
        idx->second.emplace(std::move(mm));
      } else {
        _block_index.emplace(key, phmap::btree_set<MatrixBlockMetadata>{std::move(mm)});
      }
    }
  } catch (const std::exception &e) {
    early_return = true;
    throw std::runtime_error("an error occurred in write_compressed_blocks thread: " +
                             std::string{e.what()});

  } catch (...) {
    early_return = true;
    throw;
  }
}

}  // namespace hictk::hic::internal
