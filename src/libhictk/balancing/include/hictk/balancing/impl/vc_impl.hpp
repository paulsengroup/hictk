// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

#include "hictk/chromosome.hpp"
#include "hictk/transformers/pixel_merger.hpp"

namespace hictk::balancing {

template <typename File>
inline VC::VC(const File& f, Type type) {
  switch (type) {
    case Type::cis: {
      auto res = compute_cis(f);
      _chrom_offsets = std::move(res.offsets);
      _biases = std::move(res.weights);
      _scale = std::move(res.scales);
      return;
    }
    case Type::trans: {
      auto res = compute_trans(f);
      _chrom_offsets = std::move(res.offsets);
      _biases = std::move(res.weights);
      _scale = std::move(res.scales);
      return;
    }
    case Type::gw: {
      auto res = compute_gw(f);
      _chrom_offsets = std::move(res.offsets);
      _biases = std::move(res.weights);
      _scale = std::move(res.scales);
    }
  }
}

template <typename PixelIt>
inline VC::VC(PixelIt first, PixelIt last, const hictk::BinTable& bins) {
  using N = decltype(first->count);

  const auto offset = bins.num_bin_prefix_sum().front();
  const auto size = bins.num_bin_prefix_sum().back() - offset;
  _biases.resize(size, 0);

  std::for_each(first, last, [&](const ThinPixel<N>& p) {
    _biases[p.bin1_id - offset] += p.count;
    if (p.bin1_id != p.bin2_id) {
      _biases[p.bin2_id - offset] += p.count;
    }
  });

  double sum = 0;
  double norm_sum = 0;

  std::for_each(first, last, [&](const ThinPixel<N>& p) {
    sum += p.count;
    norm_sum += p.count / (_biases[p.bin1_id - offset] * _biases[p.bin2_id - offset]);
    if (p.bin1_id != p.bin2_id) {
      sum += p.count;
      norm_sum += p.count / (_biases[p.bin1_id - offset] * _biases[p.bin2_id - offset]);
    }
  });

  _chrom_offsets.push_back(0);
  _chrom_offsets.push_back(_biases.size());
  _scale.push_back(std::sqrt(norm_sum / sum));
}

inline std::vector<double> VC::get_weights(bool rescale) const {
  if (!rescale) {
    return _biases;
  }

  std::vector<double> biases(_biases.size());
  std::uint64_t chrom_id = 0;
  for (std::size_t i = 0; i < _biases.size(); ++i) {
    if (i >= _chrom_offsets[chrom_id + 1]) {
      chrom_id++;
    }
    biases[i] = _biases[i] * _scale[chrom_id];
  }

  std::transform(biases.begin(), biases.end(), biases.begin(), [](const double n) {
    if (std::isnan(n)) {
      return 1.0;
    }
    return n;
  });

  return biases;
}

inline const std::vector<double>& VC::get_scale() const noexcept { return _scale; }

template <typename File>
inline auto VC::compute_cis(const File& f) -> Result {
  std::vector<std::uint64_t> offsets{};
  std::vector<double> scales{};
  std::vector<double> weights{};
  for (const Chromosome& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto sel = f.fetch(chrom.name());
    const VC vc{sel.template begin<double>(), sel.template end<double>(), f.bins().subset(chrom)};

    offsets.push_back(f.bins().subset(chrom).num_bin_prefix_sum().front());

    const auto chrom_weights = vc.get_weights(false);
    scales.push_back(vc.get_scale().front());
    weights.insert(weights.end(), chrom_weights.begin(), chrom_weights.end());
  }

  offsets.push_back(f.bins().size());

  return {offsets, scales, weights};
}

template <typename File>
inline auto VC::compute_trans(const File& f) -> Result {
  using PixelIt = decltype(f.fetch("chr1").template begin<double>());
  std::vector<PixelIt> heads{};
  std::vector<PixelIt> tails{};
  for (const Chromosome& chrom1 : f.chromosomes()) {
    if (chrom1.is_all()) {
      continue;
    }
    for (auto chrom2_id = chrom1.id() + 1; chrom2_id < f.chromosomes().size(); ++chrom2_id) {
      const auto& chrom2 = f.chromosomes().at(chrom2_id);
      const auto sel = f.fetch(chrom1.name(), chrom2.name());
      heads.emplace_back(sel.template begin<double>());
      tails.emplace_back(sel.template end<double>());
    }
  }

  const auto sel = transformers::PixelMerger(heads, tails);
  const VC vc{sel.begin(), sel.end(), f.bins()};

  return {{0, f.bins().size()}, vc.get_scale(), vc.get_weights(false)};
}

template <typename File>
inline auto VC::compute_gw(const File& f) -> Result {
  const auto sel = f.fetch();
  const VC vc{sel.template begin<double>(), sel.template end<double>(), f.bins()};

  return {{0, f.bins().size()}, vc.get_scale(), vc.get_weights(false)};
}

}  // namespace hictk::balancing
