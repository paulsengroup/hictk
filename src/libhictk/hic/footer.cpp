// Copyright (C) 2026 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/hic/footer.hpp"

#include <cassert>
#include <cstddef>
#include <functional>
#include <memory>
#include <utility>

#include "hictk/balancing/methods.hpp"
#include "hictk/hash.hpp"
#include "hictk/hic/index.hpp"
#include "hictk/weights.hpp"

namespace hictk::hic::internal {

bool HiCFooterMetadata::operator==(const HiCFooterMetadata &other) const noexcept {
  return url == other.url && matrix_type == other.matrix_type &&
         normalization == other.normalization && unit == other.unit &&
         resolution == other.resolution && chrom1 == other.chrom1 && chrom2 == other.chrom2;
}

bool HiCFooterMetadata::operator!=(const HiCFooterMetadata &other) const noexcept {
  return !(*this == other);
}

HiCFooter::HiCFooter(Index index_, HiCFooterMetadata metadata_, std::vector<double> expected_values,
                     std::shared_ptr<Weights> weights1, std::shared_ptr<Weights> weights2) noexcept
    : _index(std::move(index_)),
      _metadata(std::move(metadata_)),
      _expectedValues(std::move(expected_values)),
      _weights1(std::move(weights1)),
      _weights2(std::move(weights2)) {}

bool HiCFooter::operator==(const HiCFooter &other) const noexcept {
  return metadata() == other.metadata();
}

bool HiCFooter::operator!=(const HiCFooter &other) const noexcept { return !(*this == other); }

const Index &HiCFooter::index() const noexcept { return _index; }

balancing::Method HiCFooter::normalization() const noexcept { return metadata().normalization; }

const Weights &HiCFooter::weights1() const noexcept {
  assert(_weights1);
  return *_weights1;
}

const Weights &HiCFooter::weights2() const noexcept {
  assert(_weights2);
  return *_weights2;
}

}  // namespace hictk::hic::internal

std::size_t std::hash<hictk::hic::internal::HiCFooterMetadata>::operator()(
    hictk::hic::internal::HiCFooterMetadata const &m) const noexcept {
  return hictk::internal::hash_combine(0, m.url, m.matrix_type, m.normalization.to_string(), m.unit,
                                       m.resolution, m.chrom1, m.chrom2);
}

std::size_t std::hash<hictk::hic::internal::HiCFooter>::operator()(
    hictk::hic::internal::HiCFooter const &f) const noexcept {
  return std::hash<hictk::hic::internal::HiCFooterMetadata>{}(f.metadata());
}
