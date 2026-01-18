// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/pixel_selector.hpp"

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/balancing/weights.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/hic/cache.hpp"
#include "hictk/hic/common.hpp"
#include "hictk/hic/file_reader.hpp"
#include "hictk/hic/footer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic {

PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                             std::shared_ptr<const internal::HiCFooter> footer_,
                             std::shared_ptr<internal::BlockCache> cache_,
                             std::shared_ptr<const BinTable> bins_, const PixelCoordinates &coords,
                             std::optional<std::uint64_t> diagonal_band_width)
    : PixelSelector(std::move(hfs_), std::move(footer_), std::move(cache_), std::move(bins_),
                    coords, coords, diagonal_band_width) {}

PixelSelector::PixelSelector(std::shared_ptr<internal::HiCFileReader> hfs_,
                             std::shared_ptr<const internal::HiCFooter> footer_,
                             std::shared_ptr<internal::BlockCache> cache_,
                             std::shared_ptr<const BinTable> bins_, PixelCoordinates coord1_,
                             PixelCoordinates coord2_,
                             std::optional<std::uint64_t> diagonal_band_width)
    : _reader(std::make_shared<internal::HiCBlockReader>(std::move(hfs_), footer_->index(),
                                                         std::move(bins_), std::move(cache_))),
      _footer(std::move(footer_)),
      _coord1(std::make_shared<const PixelCoordinates>(std::move(coord1_))),
      _coord2(std::make_shared<const PixelCoordinates>(std::move(coord2_))),
      _diagonal_band_width(diagonal_band_width) {
  const auto query_is_cis = _coord1->bin1.chrom() == _coord2->bin1.chrom();
  if ((!query_is_cis && _coord1->bin1 > _coord2->bin1) ||
      (query_is_cis && _coord1->bin1.start() > _coord2->bin1.start())) {
    throw std::runtime_error(fmt::format(
        FMT_STRING("query {}:{}-{}; {}:{}-{}; overlaps with the lower-triangle of the matrix"),
        _coord1->bin1.chrom().name(), _coord1->bin1.start(), _coord1->bin2.end(),
        _coord2->bin1.chrom().name(), _coord2->bin1.start(), _coord2->bin2.end()));
  }
}

// NOLINTNEXTLINE(bugprone-exception-escape)
PixelSelector::~PixelSelector() noexcept {
  try {
    if (_reader) {
      clear_cache();
    }
  } catch (const std::exception &e) {
    SPDLOG_WARN(
        FMT_STRING("failed to clear the PixelSelector interaction cache of file \"{}\". "
                   "Applications are not affected. However, if this happens repeatedly the system "
                   "may run out of available memory. Reason: {}"),
        _footer->path(), e.what());
  } catch (...) {
    SPDLOG_WARN(
        FMT_STRING("failed to clear the PixelSelector interaction cache of file \"{}\". "
                   "Applications are not affected. However, if this happens repeatedly the system "
                   "may run out of available memory."),
        _footer->path());
  }
}

bool PixelSelector::operator==(const PixelSelector &other) const noexcept {
  return _reader->index().chrom1() == _reader->index().chrom2() && *_coord1 == *other._coord1 &&
         *_coord2 == *other._coord2;
}

bool PixelSelector::operator!=(const PixelSelector &other) const noexcept {
  return !(*this == other);
}

PixelSelector::PixelSelector(std::shared_ptr<internal::HiCBlockReader> reader_,
                             std::shared_ptr<const internal::HiCFooter> footer_,
                             std::shared_ptr<const PixelCoordinates> coord1_,
                             std::shared_ptr<const PixelCoordinates> coord2_,
                             std::optional<std::uint64_t> diagonal_band_width)
    : _reader(std::move(reader_)),
      _footer(std::move(footer_)),
      _coord1(std::move(coord1_)),
      _coord2(std::move(coord2_)),
      _diagonal_band_width(diagonal_band_width) {}

PixelSelector PixelSelector::fetch(PixelCoordinates coord1_, PixelCoordinates coord2_) const {
  return {_reader, _footer, std::make_shared<const PixelCoordinates>(std::move(coord1_)),
          std::make_shared<const PixelCoordinates>(std::move(coord2_)), _diagonal_band_width};
}

const PixelCoordinates &PixelSelector::coord1() const noexcept { return *_coord1; }
const PixelCoordinates &PixelSelector::coord2() const noexcept { return *_coord2; }
std::uint64_t PixelSelector::size(bool upper_triangle) const {
  if (!_coord1) {
    assert(!_coord2);
    return 0;
  }

  assert(bins().type() == BinTable::Type::fixed);

  if (_coord1->bin1.chrom() != _coord2->bin1.chrom()) {
    const auto height = _coord1->bin2.id() - _coord1->bin1.id() + 1;
    const auto width = _coord2->bin2.id() - _coord2->bin1.id() + 1;
    return height * width;
  }

  const auto start1 = _coord1->bin1.start();
  const auto start2 = _coord2->bin1.start();

  const auto end1 = ((_coord1->bin2.rel_id() + 1) * resolution()) + 1;
  const auto end2 = ((_coord2->bin2.rel_id() + 1) * resolution()) + 1;

  return area(start1, end1, start2, end2, resolution(), upper_triangle);
}

MatrixType PixelSelector::matrix_type() const noexcept { return metadata().matrix_type; }
const balancing::Method &PixelSelector::normalization() const noexcept {
  return metadata().normalization;
}
MatrixUnit PixelSelector::unit() const noexcept { return _reader->index().unit(); }
std::uint32_t PixelSelector::resolution() const noexcept {
  assert(_footer);
  return _footer->resolution();
}

const Chromosome &PixelSelector::chrom1() const noexcept { return _coord1->bin1.chrom(); }
const Chromosome &PixelSelector::chrom2() const noexcept { return _coord2->bin1.chrom(); }

const balancing::Weights &PixelSelector::weights1() const noexcept { return _footer->weights1(); }
const balancing::Weights &PixelSelector::weights2() const noexcept { return _footer->weights2(); }

const BinTable &PixelSelector::bins() const noexcept { return _reader->bins(); }
std::shared_ptr<const BinTable> PixelSelector::bins_ptr() const noexcept {
  return _reader->bins_ptr();
}

const internal::HiCFooterMetadata &PixelSelector::metadata() const noexcept {
  assert(!!_footer);
  return _footer->metadata();
}

bool PixelSelector::is_inter() const noexcept { return !is_intra(); }

bool PixelSelector::is_intra() const noexcept { return chrom1() == chrom2(); }

bool PixelSelector::empty() const noexcept { return _reader->index().empty(); }

std::size_t PixelSelector::estimate_optimal_cache_size(
    [[maybe_unused]] std::size_t num_samples) const {
  if (_reader->index().empty()) {
    return 0;  // should we throw instead?
  }

  std::seed_seq sseq({_reader->index().size()});
  std::mt19937_64 rand_eng(sseq);

  // Try to guess the average block size
  std::size_t avg_block_size = 0;
  const auto &idx = _reader->index();

  std::vector<internal::BlockIndex> blocks(std::min(idx.size(), num_samples));
  std::sample(idx.begin(), idx.end(), blocks.begin(), blocks.size(), rand_eng);
  for (const auto &blki : blocks) {
    avg_block_size += _reader->read_size(chrom1(), chrom2(), blki);
  }
  avg_block_size /= blocks.size();

  // Try to guess how many blocks overlap a single row of pixels
  std::size_t max_blocks_per_row = 0;
  const auto &chrom = coord1().bin1.chrom();
  const auto bin_size = bins().resolution();

  const std::size_t first_bin_id = 0;
  const std::size_t last_bin_id = bins().at(chrom, chrom.size() - 1).rel_id();

  const auto samples = std::min(num_samples, bins().subset(chrom).size());
  for (std::size_t i = 0; i < samples; ++i) {
    const auto bin_id = std::uniform_int_distribution<std::size_t>{
        first_bin_id, std::min(last_bin_id, last_bin_id - 1)}(rand_eng);

    const auto pos1 = static_cast<std::uint32_t>(bin_id * bin_size);
    const auto bin1 = bins().at(chrom, pos1);

    auto overlap = idx.find_overlaps({bin1, bin1}, coord2(), _diagonal_band_width);
    const auto num_blocks = static_cast<std::size_t>(std::distance(overlap.begin(), overlap.end()));
    max_blocks_per_row = (std::max)(max_blocks_per_row, num_blocks);
  }

  return max_blocks_per_row * avg_block_size * sizeof(ThinPixel<float>);
}

void PixelSelector::clear_cache() const {
  if (_reader->cache_size() == 0) {
    return;
  }
  for (auto blki : _reader->index().find_overlaps(coord1(), coord2(), _diagonal_band_width)) {
    _reader->evict(coord1().bin1.chrom(), coord2().bin1.chrom(), blki);
  }
}

PixelSelectorAll::PixelSelectorAll(std::vector<PixelSelector> selectors_,
                                   std::shared_ptr<internal::WeightCache> weight_cache)
    : _selectors(std::move(selectors_)),
      _bins(_selectors.empty() ? nullptr : _selectors.front().bins_ptr()),
      _weight_cache(std::move(weight_cache)) {
  if (_selectors.empty()) {
    throw std::invalid_argument(
        "selectors_ cannot be empty! You may want to call a different constructor.");
  }

  assert(!!_bins);
}

PixelSelectorAll::PixelSelectorAll(std::shared_ptr<const BinTable> bins_,
                                   std::shared_ptr<internal::WeightCache> weight_cache) noexcept
    : _bins(std::move(bins_)), _weight_cache(std::move(weight_cache)) {}

bool PixelSelectorAll::empty() const noexcept { return begin<float>() == end<float>(); }

std::uint64_t PixelSelectorAll::size(bool upper_triangle) const noexcept {
  const auto n = bins().size();
  if (upper_triangle) {
    return (n * (n + 1)) / 2;
  }
  return n * n;
}

MatrixType PixelSelectorAll::matrix_type() const noexcept {
  return _selectors.front().matrix_type();
}
const balancing::Method &PixelSelectorAll::normalization() const noexcept {
  return _selectors.front().normalization();
}
MatrixUnit PixelSelectorAll::unit() const noexcept { return _selectors.front().unit(); }
std::uint32_t PixelSelectorAll::resolution() const noexcept { return bins().resolution(); }

const BinTable &PixelSelectorAll::bins() const noexcept {
  assert(_bins);
  return *_bins;
}
std::shared_ptr<const BinTable> PixelSelectorAll::bins_ptr() const noexcept { return _bins; }

const balancing::Weights &PixelSelectorAll::weights() const {
  if (!_weight_cache) {
    throw std::runtime_error(
        "PixelSelectorAll::weights() was called on an instance of PixelSelectorAll with null "
        "_weight_cache");
  }
  auto weights = _weight_cache->get_or_init(0, normalization());
  if (!weights->empty()) {
    return *weights;
  }

  if (normalization() == balancing::Method::NONE()) {
    *weights = balancing::Weights{1.0, bins().size(), balancing::Weights::Type::DIVISIVE};
    return *weights;
  }

  std::vector<double> buff(bins().size(), std::numeric_limits<double>::quiet_NaN());
  std::for_each(_selectors.begin(), _selectors.end(), [&](const PixelSelector &sel) {
    if (sel.is_intra()) {
      const auto &chrom_weights = sel.weights1();
      const auto offset = static_cast<std::ptrdiff_t>(bins().at(sel.chrom1()).id());
      std::copy(chrom_weights.begin(balancing::Weights::Type::DIVISIVE),
                chrom_weights.end(balancing::Weights::Type::DIVISIVE), buff.begin() + offset);
    }
  });

  *weights = balancing::Weights{std::move(buff), balancing::Weights::Type::DIVISIVE};
  return *weights;
}

}  // namespace hictk::hic
