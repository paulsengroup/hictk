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

inline std::string MatrixMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(chr1Idx);
  buffer.write(chr2Idx);
  buffer.write(nResolutions);

  return buffer.get();
}

inline std::string MatrixBlockMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(blockNumber);
  buffer.write(blockPosition);
  buffer.write(blockSizeBytes);

  return buffer.get();
}

inline bool MatrixBlockMetadata::operator<(const MatrixBlockMetadata &other) const noexcept {
  return blockNumber < other.blockNumber;
}

inline std::string MatrixResolutionMetadata::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

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
                                                     std::string &compression_buffer,
                                                     bool clear) const {
  // TODO support dense layout
  // TODO support representation using shorts

  if (clear) {
    buffer.clear();
  }

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

inline std::string MasterIndex::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(key);
  buffer.write(position);
  buffer.write(size);

  return buffer.get();
}

inline ExpectedValuesBlock::ExpectedValuesBlock(std::string_view unit_, std::uint32_t bin_size,
                                                const std::vector<double> &weights,
                                                const std::vector<std::uint32_t> &chrom_ids,
                                                const std::vector<double> &scale_factors)
    : unit(std::string{unit_}),
      binSize(static_cast<std::int32_t>(bin_size)),
      nValues(static_cast<std::int32_t>(weights.size())),
      value(weights.size()),
      nChrScaleFactors(static_cast<std::int32_t>(chrom_ids.size())),
      chrIndex(chrom_ids.size()),
      chrScaleFactor(chrom_ids.size()) {
  std::transform(weights.begin(), weights.end(), value.begin(),
                 [](const auto n) { return static_cast<float>(n); });
  std::transform(chrom_ids.begin(), chrom_ids.end(), chrIndex.begin(),
                 [](const auto n) { return static_cast<std::int32_t>(n); });
  std::transform(scale_factors.begin(), scale_factors.end(), chrScaleFactor.begin(),
                 [](const auto n) { return static_cast<float>(n); });
}

inline std::string ExpectedValuesBlock::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(unit);
  buffer.write(binSize);
  buffer.write(nValues);
  buffer.write(value);
  buffer.write(nChrScaleFactors);

  assert(chrIndex.size() == chrScaleFactor.size());
  for (std::size_t i = 0; i < chrIndex.size(); ++i) {
    buffer.write(chrIndex[i]);
    buffer.write(chrScaleFactor[i]);
  }

  return buffer.get();
}

inline std::string ExpectedValues::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(nExpectedValueVectors);

  if (nExpectedValueVectors == 0) {
    return buffer.get();
  }

  for (const auto &ev : expectedValues) {
    std::ignore = ev.serialize(buffer, false);
  }

  return buffer.get();
}

inline std::string NormalizedExpectedValues::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(nNormExpectedValueVectors);

  return buffer.get();
}

inline std::string NormalizationVectorIndex::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(nNormVectors);

  return buffer.get();
}

inline std::string NormalizationVectorArray::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

  buffer.write(nValues);

  return buffer.get();
}

inline std::string FooterV5::serialize(BinaryBuffer &buffer, bool clear) const {
  if (clear) {
    buffer.clear();
  }

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

inline ExpectedValuesAggregator::ExpectedValuesAggregator(std::shared_ptr<const BinTable> bins)
    : _bins(std::move(bins)) {
  std::uint32_t max_length = 0;
  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    if (chrom1.is_all()) {
      continue;
    }

    max_length = std::max(max_length, chrom1.size());
    _cis_sum.emplace(chrom1, 0.0);

    for (std::uint32_t chrom2_id = chrom1_id + 1; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);
      _trans_sum.emplace(std::make_pair(chrom1, chrom2), 0.0);
    }

    for (const auto &chrom : chromosomes()) {
      if (chrom.is_all()) {
        continue;
      }
      _num_bins_gw += chrom.size();
    }
  }

  const auto bin_size = _bins->bin_size();
  const auto max_n_bins = (max_length + bin_size - 1) / bin_size;
  _possible_distances.resize(max_n_bins, 0.0);
  _actual_distances.resize(max_n_bins, 0.0);

  for (const auto &chrom : chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto n_bins = chrom.size() / bin_size;
    for (std::uint32_t i = 0; i < n_bins; ++i) {
      _possible_distances[i] += n_bins - i;
    }
  }
}

inline void ExpectedValuesAggregator::add(const ThinPixel<float> &p) {
  add(Pixel<float>{*_bins, p});
}

inline void ExpectedValuesAggregator::add(const Pixel<float> &p) {
  const auto &chrom1 = p.coords.bin1.chrom();
  const auto &chrom2 = p.coords.bin2.chrom();

  if (p.coords.is_intra()) {
    at(chrom1) += static_cast<double>(p.count);
    const auto i = p.coords.bin2.id() - p.coords.bin1.id();
    _actual_distances[i] += static_cast<double>(p.count);
  } else {
    at(chrom1, chrom2) += static_cast<double>(p.count);
  }
}

inline void ExpectedValuesAggregator::compute_density() {
  compute_density_cis();
  compute_density_trans();
}

inline const std::vector<double> &ExpectedValuesAggregator::weights() const noexcept {
  return _weights;
}

inline double ExpectedValuesAggregator::scaling_factor(const Chromosome &chrom) const {
  return _scaling_factors.at(chrom);
}

inline const phmap::btree_map<Chromosome, double> &ExpectedValuesAggregator::scaling_factors()
    const noexcept {
  return _scaling_factors;
}

inline void ExpectedValuesAggregator::compute_density_cis() {
  // Re-implementation of the algorithm used by HiCTools:
  // https://github.com/aidenlab/HiCTools/blob/6b2fab8e78685deae199c33bbb167dcab1dbfbb3/src/hic/tools/utils/original/ExpectedValueCalculation.java#L184

  auto num_sum = _actual_distances.front();
  auto den_sum = _possible_distances.front();
  std::size_t bound1 = 0;
  std::size_t bound2 = 0;

  const auto shot_noise_minimum = 400.0;
  const auto max_num_bins = _actual_distances.size();

  _weights.resize(max_num_bins);
  std::fill(_weights.begin(), _weights.end(), 0.0);

  for (std::size_t ii = 0; ii < max_num_bins; ii++) {
    if (num_sum < shot_noise_minimum) {
      while (num_sum < shot_noise_minimum && bound2 < max_num_bins) {
        bound2++;
        num_sum += _actual_distances[bound2];
        den_sum += _possible_distances[bound2];
      }
    } else if (num_sum >= shot_noise_minimum && bound2 - bound1 > 0) {
      while (bound2 > bound1 && bound2 < _num_bins_gw && bound1 < _num_bins_gw &&
             num_sum - _actual_distances[bound1] - _actual_distances[bound2] >=
                 shot_noise_minimum) {
        num_sum = num_sum - _actual_distances[bound1] - _actual_distances[bound2];
        den_sum = den_sum - _possible_distances[bound1] - _possible_distances[bound2];
        bound1++;
        bound2--;
      }
    }
    _weights[ii] = num_sum / den_sum;
    if (bound2 + 2 < max_num_bins) {
      num_sum += _actual_distances[bound2 + 1] + _actual_distances[bound2 + 2];
      den_sum += _possible_distances[bound2 + 1] + _possible_distances[bound2 + 2];
      bound2 += 2;
    } else if (bound2 + 1 < max_num_bins) {
      num_sum += _actual_distances[bound2 + 1];
      den_sum += _possible_distances[bound2 + 1];
      bound2++;
    }
  }

  for (const auto &chrom : chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    auto num_chrom_bins = chrom.size() / _bins->bin_size();
    auto expected_count = 0.0;
    for (std::size_t n = 0; n < num_chrom_bins; n++) {
      if (n < max_num_bins) {
        double v = _weights[n];
        expected_count += (double(num_chrom_bins) - double(n)) * v;
      }
    }

    double observed_count = _cis_sum.at(chrom);
    double f = expected_count / observed_count;
    _scaling_factors.emplace(chrom, f);
  }
}

inline void ExpectedValuesAggregator::compute_density_trans() {
  for (auto &[k, v] : _trans_sum) {
    // We round-down to match HiCTools behavior
    const auto num_bins1 = k.first.size() / _bins->bin_size();
    const auto num_bins2 = k.second.size() / _bins->bin_size();
    const auto num_pixels = num_bins1 * num_bins2;
    v = num_pixels != 0 ? v / static_cast<double>(num_pixels) : 0.0;
  }
}

inline double ExpectedValuesAggregator::at(const Chromosome &chrom) const {
  return _cis_sum.at(chrom);
}

inline double ExpectedValuesAggregator::at(const Chromosome &chrom1,
                                           const Chromosome &chrom2) const {
  return _trans_sum.at(std::make_pair(chrom1, chrom2));
}

inline double &ExpectedValuesAggregator::at(const Chromosome &chrom) { return _cis_sum.at(chrom); }

inline double &ExpectedValuesAggregator::at(const Chromosome &chrom1, const Chromosome &chrom2) {
  return _trans_sum.at(std::make_pair(chrom1, chrom2));
}

inline const Reference &ExpectedValuesAggregator::chromosomes() const noexcept {
  assert(_bins);
  return _bins->chromosomes();
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
  for (auto &[k, v] : _pixel_tank) {
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
    auto [it, _] = _bin_tables.emplace(resolution,
                                       std::make_shared<const BinTable>(chromosomes(), resolution));
    _pixel_tank.emplace(resolution, it->second);
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
  return *_bin_tables.at(resolution);
}

inline const std::vector<std::uint32_t> &HiCFileWriter::resolutions() const noexcept {
  assert(_header);
  return _header->resolutions;
}

template <typename PixelIt, typename>
inline void HiCFileWriter::append_pixels(std::uint32_t resolution, PixelIt first_pixel,
                                         PixelIt last_pixel, bool update_expected_values) {
  _pixel_tank.at(resolution).add_pixels(first_pixel, last_pixel, update_expected_values);
}

inline void HiCFileWriter::write_pixels([[maybe_unused]] bool write_chromosome_ALL) {
  coarsen_pixels();

  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);
      for (const auto &res : resolutions()) {
        write_pixels(chrom1, chrom2, res);
      }
    }
  }
  if (write_chromosome_ALL) {
    // write_pixels_ALL();
  }
}

inline void HiCFileWriter::write_pixels(const Chromosome &chrom1, const Chromosome &chrom2,
                                        std::uint32_t resolution) {
  struct PixelBlock {
    std::vector<ThinPixel<float>> pixels;
    std::uint32_t bin1_id_offset;
    std::uint32_t bin2_id_offset;
  };

  phmap::btree_map<std::uint64_t, PixelBlock> blocks;

  if (!_pixel_tank.at(resolution).contains(chrom1, chrom2)) {
    return;
  }

  {
    const auto &bin_table = bins(resolution);
    assert(!_pixel_tank.empty());
    const auto &pixels = _pixel_tank.at(resolution).pixels(chrom1, chrom2);
    const auto num_bins = compute_num_bins(chrom1.id(), chrom2.id(), bin_table.bin_size());
    const auto num_columns =
        compute_block_column_count(num_bins, bin_table.bin_size(),
                                   chrom1 == chrom2 ? DEFAULT_INTRA_CUTOFF : DEFAULT_INTER_CUTOFF);
    const auto num_rows = num_bins / num_columns + 1;

    const BlockMapperIntra mapper_intra(num_rows, num_columns);
    const BlockMapperInter mapper_inter(num_rows, num_columns);

    for (const auto &p : pixels) {
      Pixel<float> pp(bin_table, p);

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
    write_interaction_block(block_id, chrom1, chrom2, resolution, blk.pixels, blk.bin1_id_offset,
                            blk.bin2_id_offset);
  }
}

inline void HiCFileWriter::write_pixels_ALL([[maybe_unused]] std::size_t num_bins) {
  /*
  const auto genome_size = chromosomes().chrom_size_prefix_sum().back();
  const auto resolution = static_cast<std::uint32_t>(genome_size / num_bins);

  _bin_tables.emplace(resolution, std::make_shared<const BinTable>(chromosomes(), resolution));

  const auto &bin_table = bins(resolutions().front());
  const auto factor = std::max(bin_table.size(), bin_table.size() / resolution);

  assert(!_pixel_tank.empty());
  std::vector<ThinPixel<float>> pixels{};
  for (const auto &[_, pixels_] : _pixel_tank) {
    const transformers::CoarsenPixels coarsener(pixels_.begin(), pixels_.end(),
                                                std::make_shared<BinTable>(bin_table), factor);
    std::copy(coarsener.begin(), coarsener.end(), std::back_inserter(pixels));
  }

  std::sort(pixels.begin(), pixels.end());
  write_interaction_block(0, chromosomes().at(0), chromosomes().at(0), resolution, pixels, 0, 0);
*/
}

inline void HiCFileWriter::write_header() {
  assert(_fs->tellg() == 0);

  assert(_header);
  assert(_header->version == 9);
  assert(!chromosomes().empty());
  assert(!resolutions().empty());

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
}

inline void HiCFileWriter::write_footer_offset(std::int64_t master_index_offset) {
  const auto foffset = _fs->tellp();
  const auto offset = sizeof("HIC") + sizeof(_header->version);
  _fs->seekp(offset);
  _fs->write(master_index_offset);
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
                                                      std::uint32_t chrom2_id, float sum_counts,
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
    m.sumCounts = sum_counts;
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

inline void HiCFileWriter::coarsen_pixels() {
  assert(!resolutions().empty());
  const auto base_resolution = resolutions().front();
  _pixel_tank.at(base_resolution).finalize();

  for (std::size_t i = 1; i < resolutions().size(); ++i) {
    const auto &resolution = resolutions()[i];
    assert(resolution % base_resolution == 0);
    const auto factor = resolution / resolutions().front();
    for (const auto &[k, v] : _pixel_tank.at(base_resolution).pixels()) {
      transformers::CoarsenPixels coarsener(v.begin(), v.end(), _bin_tables.at(base_resolution),
                                            factor);
      _pixel_tank.at(resolution).add_pixels(coarsener.begin(), coarsener.end());
    }
    _pixel_tank.at(resolution).finalize();
  }
}

inline auto HiCFileWriter::write_body_metadata(std::uint32_t resolution, const Chromosome &chrom1,
                                               const Chromosome &chrom2, const std::string &unit) {
  struct Result {
    std::size_t matrix_metadata_offset;
    std::size_t matrix_metadata_bytes;
  };

  const auto tot_counts = _pixel_tank.at(resolution).contains(chrom1, chrom2)
                              ? _pixel_tank.at(resolution).matrix_counts(chrom1, chrom2)
                              : 0.0F;

  Result offsets{};
  offsets.matrix_metadata_offset = write_matrix_metadata(chrom1.id(), chrom2.id());
  write_resolutions_metadata(chrom1.id(), chrom2.id(), tot_counts, unit);

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
  const std::int64_t nBytesV5 = -1;
  const std::int32_t nEntries = static_cast<std::int32_t>(_footers.size());
  _fs->write(nBytesV5);
  _fs->write(nEntries);

  for (const auto &f : _footers) {
    _fs->write(f.masterIndex.serialize(_bbuffer));
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

inline void HiCFileWriter::write_footer_section_size(std::streamoff footer_offset,
                                                     std::uint64_t bytes) {
  const auto offset = _fs->tellp();
  _fs->seekp(static_cast<std::int64_t>(footer_offset));
  _fs->write(static_cast<std::int64_t>(bytes));
  _fs->seekp(static_cast<std::int64_t>(offset));
}

inline void HiCFileWriter::write_expected_values(std::string_view unit) {
  ExpectedValues evs{};
  evs.nExpectedValueVectors = static_cast<std::int32_t>(_pixel_tank.size());

  for (const auto &[resolution, pt] : _pixel_tank) {
    std::vector<std::uint32_t> chrom_ids{};
    std::vector<double> scale_factors{};
    const auto &ev = _pixel_tank.at(resolution).expected_values();
    assert(!ev.weights().empty());
    for (const auto &[chrom, factor] : ev.scaling_factors()) {
      chrom_ids.push_back(chrom.id());
      scale_factors.push_back(factor);
    }
    evs.expectedValues.emplace_back(unit, resolution, ev.weights(), chrom_ids, scale_factors);
  }

  _fs->write(evs.serialize(_bbuffer));
}

inline void HiCFileWriter::finalize() {
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, std::size_t> file_offsets{};
  phmap::flat_hash_map<std::pair<Chromosome, Chromosome>, std::size_t> matrix_size_bytes{};

  std::size_t master_index_offset = std::numeric_limits<std::size_t>::max();

  for (std::uint32_t chrom1_id = 0; chrom1_id < chromosomes().size(); ++chrom1_id) {
    const auto &chrom1 = chromosomes().at(chrom1_id);
    for (std::uint32_t chrom2_id = chrom1_id; chrom2_id < chromosomes().size(); ++chrom2_id) {
      const auto &chrom2 = chromosomes().at(chrom2_id);  // TODO handle ALL
      for (std::size_t i = 0; i < resolutions().size(); ++i) {
        const auto res = resolutions()[i];
        const auto offsets = write_body_metadata(res, chrom1, chrom2);
        const auto key = std::make_pair(chrom1, chrom2);

        if (i == 0) {
          master_index_offset = std::min(offsets.matrix_metadata_offset, master_index_offset);
          file_offsets[key] = offsets.matrix_metadata_offset;
          matrix_size_bytes[key] = offsets.matrix_metadata_bytes;
        }
        if (chrom2.is_all()) {
          break;
        }
      }
      if (chrom2.is_all()) {
        break;
      }
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

  write_expected_values("BP");

  write_footer_section_size(footer_offset, _fs->tellp() - static_cast<std::size_t>(footer_offset));

  const auto normVectorIndexPosition = _fs->tellp();
  _fs->write(std::int32_t(0));  // no nNormExpectedValueVectors
  _fs->write(std::int32_t(0));  // no NormVectors
  const auto normVectorIndexLength = _fs->tellp() - normVectorIndexPosition;

  write_footer_offset(footer_offset);
  write_norm_vector_index(static_cast<std::streamoff>(normVectorIndexPosition),
                          normVectorIndexLength);
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
  return std::clamp(num_columns, std::size_t(1), max_sqrt - 1);
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
