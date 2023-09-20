// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cmath>
#include <iostream>
#include <iterator>
#include <nonstd/span.hpp>
#include <numeric>
#include <type_traits>

#include "hictk/cooler/cooler.hpp"
#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::balancing {

inline bool SparseMatrix::empty() const noexcept { return size() == 0; }
inline std::size_t SparseMatrix::size() const noexcept { return counts.size(); }

inline void SparseMatrix::clear() noexcept {
  bin1_ids.clear();
  bin2_ids.clear();
  counts.clear();
  chrom_offsets.clear();
}

template <typename File>
inline ICE::ICE(const File& f, Type type, double tol, std::size_t max_iters,
                std::size_t num_masked_diags, std::size_t min_nnz, std::size_t min_count,
                double mad_max)
    : _chrom_bin_offsets(read_chrom_bin_offsets(f.bins())), _biases(f.bins().size(), 1.0) {
  auto matrix = construct_sparse_matrix(f, type, num_masked_diags);

  initialize_biases(matrix.bin1_ids, matrix.bin2_ids, matrix.counts, _biases, _chrom_bin_offsets,
                    min_nnz, min_count, mad_max);

  std::vector<double> margs(_biases.size());
  if (type != Type::cis) {
    _variance.resize(1, 0);
    _scale.resize(1, std::numeric_limits<double>::quiet_NaN());

    auto weights = type == Type::trans
                       ? compute_weights_from_chromosome_sizes(f.bins(), _chrom_bin_offsets)
                       : std::vector<double>{};
    if (type == Type::trans) {
      mask_cis_interactions(matrix.bin1_ids, matrix.bin2_ids, matrix.counts, _chrom_bin_offsets);
    }

    for (std::size_t i = 0; i < max_iters; ++i) {
      const auto res =
          inner_loop(matrix.bin1_ids, matrix.bin2_ids, matrix.counts, _biases, margs, 0, weights);
      _variance[0] = res.variance;
      _scale[0] = res.scale;
      if (res.variance < tol) {
        return;
      }
    }
  }

  _variance.resize(_chrom_bin_offsets.size() - 1, 0);
  _scale.resize(_chrom_bin_offsets.size() - 1, std::numeric_limits<double>::quiet_NaN());
  for (std::size_t i = 1; i < _chrom_bin_offsets.size(); ++i) {
    const auto i0 = matrix.chrom_offsets[i - 1];
    const auto i1 = matrix.chrom_offsets[i];

    auto bin1_ids_ = nonstd::span(matrix.bin1_ids).subspan(i0, i1 - i0);
    auto bin2_ids_ = nonstd::span(matrix.bin2_ids).subspan(i0, i1 - i0);
    std::vector<double> counts_{};
    std::copy(matrix.counts.begin() + static_cast<std::ptrdiff_t>(i0),
              matrix.counts.begin() + static_cast<std::ptrdiff_t>(i1), std::back_inserter(counts_));

    const auto j0 = _chrom_bin_offsets[i - 1];
    const auto j1 = _chrom_bin_offsets[i];

    auto biases_ = nonstd::span(_biases).subspan(j0, j1 - j0);
    auto margs_ = nonstd::span(margs).subspan(j0, j1 - j0);
    for (std::size_t k = 0; k < max_iters; ++k) {
      const auto res = inner_loop(bin1_ids_, bin2_ids_, counts_, biases_, margs_, j0);
      _variance[i - 1] = res.variance;
      _scale[i - 1] = res.scale;

      if (res.variance < tol) {
        break;
      }
    }
  }
}

template <typename File>
auto ICE::construct_sparse_matrix(const File& f, Type type, std::size_t num_masked_diags,
                                  std::size_t bin_offset) -> SparseMatrix {
  switch (type) {
    case Type::cis:
      return construct_sparse_matrix_cis(f, num_masked_diags, bin_offset);
    case Type::trans:
      [[fallthrough]];
    case Type::gw:
      return construct_sparse_matrix_gw(f, num_masked_diags, bin_offset);
  }
}

template <typename File>
inline auto ICE::construct_sparse_matrix_gw(const File& f, std::size_t num_masked_diags,
                                            std::size_t bin_offset) -> SparseMatrix {
  SparseMatrix m{};

  const auto sel = f.fetch();
  std::for_each(sel.template begin<double>(), sel.template end<double>(), [&](const auto& p) {
    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
      m.bin1_ids.push_back(p.bin1_id - bin_offset);
      m.bin2_ids.push_back(p.bin2_id - bin_offset);
      m.counts.push_back(p.count);
    }
  });

  m.bin1_ids.shrink_to_fit();
  m.bin2_ids.shrink_to_fit();
  m.counts.shrink_to_fit();

  return m;
}

template <typename File>
[[nodiscard]] inline auto ICE::construct_sparse_matrix_cis(const File& f,
                                                           std::size_t num_masked_diags,
                                                           std::size_t bin_offset) -> SparseMatrix {
  SparseMatrix m{};
  m.chrom_offsets.push_back(0);

  for (const Chromosome& chrom : f.chromosomes()) {
    const auto sel = f.fetch(chrom.name());
    std::for_each(sel.template begin<double>(), sel.template end<double>(),
                  [&](const ThinPixel<double>& p) {
                    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
                      m.bin1_ids.push_back(p.bin1_id - bin_offset);
                      m.bin2_ids.push_back(p.bin2_id - bin_offset);

                      m.counts.push_back(p.count);
                    }
                  });

    m.chrom_offsets.push_back(m.bin1_ids.size());
  }
  m.bin1_ids.shrink_to_fit();
  m.bin2_ids.shrink_to_fit();
  m.counts.shrink_to_fit();

  return m;
}

inline void ICE::marginalize(nonstd::span<const std::size_t> bin1_ids,
                             nonstd::span<const std::size_t> bin2_ids,
                             nonstd::span<const double> counts, nonstd::span<double> marg,
                             std::size_t bin_offset) {
  std::fill(marg.begin(), marg.end(), 0);

  for (std::size_t i = 0; i < counts.size(); ++i) {
    const auto i1 = bin1_ids[i] - bin_offset;
    const auto i2 = bin2_ids[i] - bin_offset;

    marg[i1] += counts[i];
    marg[i2] += counts[i];
  }
}

inline void ICE::times_outer_product_marg(nonstd::span<const std::size_t> bin1_ids,
                                          nonstd::span<const std::size_t> bin2_ids,
                                          nonstd::span<const double> counts,
                                          nonstd::span<const double> biases,
                                          nonstd::span<double> marg, std::size_t bin_offset,
                                          nonstd::span<const double> weights) {
  assert(bin1_ids.size() == counts.size());
  assert(bin2_ids.size() == counts.size());

  assert(biases.size() == marg.size());
  assert(biases.size() == weights.size() || weights.empty());

  std::fill(marg.begin(), marg.end(), 0);

  for (std::size_t i = 0; i < counts.size(); ++i) {
    const auto i1 = bin1_ids[i] - bin_offset;
    const auto i2 = bin2_ids[i] - bin_offset;
    const auto w1 = weights.empty() ? 1 : weights[i1];
    const auto w2 = weights.empty() ? 1 : weights[i2];
    const auto count = counts[i] * (w1 * biases[i1]) * (w2 * biases[i2]);

    marg[i1] += count;
    marg[i2] += count;
  }
}

inline void ICE::marginalize_nnz(nonstd::span<const std::size_t> bin1_ids,
                                 nonstd::span<const std::size_t> bin2_ids,
                                 nonstd::span<const double> counts, nonstd::span<double> marg,
                                 std::size_t bin_offset) {
  std::fill(marg.begin(), marg.end(), 0);

  for (std::size_t i = 0; i < counts.size(); ++i) {
    const auto i1 = bin1_ids[i] - bin_offset;
    const auto i2 = bin2_ids[i] - bin_offset;

    marg[i1] += counts[i] != 0;
    marg[i2] += counts[i] != 0;
  }
}

inline void ICE::min_nnz_filtering(nonstd::span<const std::size_t> bin1_ids,
                                   nonstd::span<const std::size_t> bin2_ids,
                                   nonstd::span<const double> counts, nonstd::span<double> biases,
                                   std::size_t min_nnz, nonstd::span<double> marg_buff,
                                   std::size_t bin_offset) {
  marginalize_nnz(bin1_ids, bin2_ids, counts, marg_buff, bin_offset);
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (marg_buff[i] < static_cast<double>(min_nnz)) {
      biases[i] = 0;
    }
  }
}

inline void ICE::min_count_filtering(nonstd::span<double> biases, std::size_t min_count,
                                     const nonstd::span<double> marg) {
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (marg[i] < static_cast<double>(min_count)) {
      biases[i] = 0;
    }
  }
}

inline void ICE::mad_max_filtering(nonstd::span<const std::size_t> chrom_offsets,
                                   nonstd::span<double> biases, std::vector<double> marg,
                                   double mad_max) {
  auto median = [](auto v) {
    assert(!v.empty());

    const auto size = static_cast<std::ptrdiff_t>(v.size());
    auto first = v.begin();
    auto mid = first + (size / 2);
    auto last = v.end();

    std::nth_element(first, mid, last);

    if (size % 2 != 0) {
      return *mid;
    }

    const auto n1 = *mid;
    std::nth_element(first, --mid, last);
    const auto n2 = *mid;

    return (n1 + n2) / 2;
  };

  auto mad = [&](const auto vin) {
    const auto median_ = median(vin);
    auto vout = vin;

    std::transform(vout.begin(), vout.end(), vout.begin(),
                   [&](const auto n) { return std::abs(n - median_); });

    return median(vout);
  };

  assert(chrom_offsets.size() > 1);
  std::vector<double> cmarg{};
  for (std::size_t i = 1; i < chrom_offsets.size(); ++i) {
    const auto i0 = static_cast<std::ptrdiff_t>(chrom_offsets[i - 1] - chrom_offsets.front());
    const auto i1 = static_cast<std::ptrdiff_t>(chrom_offsets[i] - chrom_offsets.front());

    cmarg.clear();
    std::copy_if(marg.begin() + i0, marg.begin() + i1, std::back_inserter(cmarg),
                 [](const auto n) { return n > 0; });

    if (!cmarg.empty()) {
      const auto median_ = median(cmarg);
      std::transform(marg.begin() + i0, marg.begin() + i1, marg.begin() + i0,
                     [&](const auto n) { return n / median_; });
    }
  }

  std::vector<double> log_nz_marg{};
  for (const auto n : marg) {
    if (n > 0) {
      log_nz_marg.push_back(std::log(n));
    }
  }

  if (log_nz_marg.empty()) {
    return;
  }

  const auto median_log_nz_marg = median(log_nz_marg);
  const auto dev_log_nz_marg = mad(log_nz_marg);

  const auto cutoff = std::exp(median_log_nz_marg - mad_max * dev_log_nz_marg);

  for (std::size_t i = 0; i < marg.size(); ++i) {
    if (marg[i] < cutoff) {
      biases[i] = 0.0;
    }
  }
}

inline auto ICE::inner_loop(nonstd::span<const std::size_t> bin1_ids,
                            nonstd::span<const std::size_t> bin2_ids,
                            nonstd::span<const double> counts, nonstd::span<double> biases,
                            nonstd::span<double> marg_buffer, std::size_t bin_offset,
                            nonstd::span<double> weights) -> Result {
  times_outer_product_marg(bin1_ids, bin2_ids, counts, biases, marg_buffer, bin_offset, weights);

  double marg_sum = 0.0;
  std::size_t nnz_marg{};
  for (const auto& n : marg_buffer) {
    marg_sum += n;
    nnz_marg += n != 0;
  }

  if (nnz_marg == 0) {
    std::fill(biases.begin(), biases.end(), std::numeric_limits<double>::quiet_NaN());
    return {std::numeric_limits<double>::quiet_NaN(), 0.0};
  }

  const auto avg_nzmarg = (marg_sum / static_cast<double>(nnz_marg));
  for (std::size_t i = 0; i < biases.size(); ++i) {
    const auto n = marg_buffer[i] / avg_nzmarg;
    if (n != 0) {
      biases[i] /= n;
    }
  }

  double ssq_nzmarg = 0.0;
  for (const auto n : marg_buffer) {
    if (n != 0) {
      ssq_nzmarg += std::pow(n - avg_nzmarg, 2);
    }
  }
  const auto var_nzmarg = ssq_nzmarg / static_cast<double>(nnz_marg - 1);

  return {avg_nzmarg, var_nzmarg};
}

inline void ICE::initialize_biases(nonstd::span<const std::size_t> bin1_ids,
                                   nonstd::span<const std::size_t> bin2_ids,
                                   nonstd::span<const double> counts, nonstd::span<double> biases,
                                   nonstd::span<const std::size_t> chrom_bin_offsets,
                                   std::size_t min_nnz, std::size_t min_count, double mad_max) {
  std::vector<double> margs(biases.size());
  if (min_nnz != 0) {
    min_nnz_filtering(bin1_ids, bin2_ids, counts, biases, min_nnz, margs);
  }

  marginalize(bin1_ids, bin2_ids, counts, margs);
  if (min_count != 0) {
    min_count_filtering(biases, min_count, margs);
  }

  if (mad_max != 0) {
    mad_max_filtering(chrom_bin_offsets, biases, margs, mad_max);
  }
}

inline std::vector<size_t> ICE::read_chrom_bin_offsets(const BinTable& bins) {
  std::vector<std::size_t> buff{0};
  for (const Chromosome& chrom : bins.chromosomes()) {
    const auto nbins = (chrom.size() + bins.bin_size() - 1) / bins.bin_size();
    buff.push_back(buff.back() + nbins);
  }

  return buff;
}

inline std::vector<double> ICE::compute_weights_from_chromosome_sizes(
    const BinTable& bins, nonstd::span<std::size_t> chrom_bin_offsets) {
  std::vector<double> weights(bins.size());
  for (std::uint32_t i = 1; i < chrom_bin_offsets.size(); ++i) {
    const auto& chrom = bins.chromosomes().at(i - 1);
    const auto i0 = chrom_bin_offsets[i - 1];
    const auto i1 = chrom_bin_offsets[i];

    const auto nbins = static_cast<double>(bins.size());
    const auto cnbins =
        std::ceil(static_cast<double>(chrom.size()) / static_cast<double>(bins.bin_size()));

    for (std::size_t j = i0; j < i1; ++j) {
      weights[j] = 1.0 / (1.0 - cnbins / nbins);
    }
  }
  return weights;
}

inline std::vector<double> ICE::get_weights([[maybe_unused]] bool rescale) const {
  std::vector<double> biases(_biases.size());
  if (!rescale) {
    return biases;
  }

  if (_scale.size() == 1) {
    const auto scale = std::sqrt(_scale[0]);
    std::transform(_biases.begin(), _biases.end(), biases.begin(), [&](const auto n) {
      return n == 0 ? std::numeric_limits<double>::quiet_NaN() : n / scale;
    });
  } else {
    for (std::size_t i = 1; i < _chrom_bin_offsets.size(); ++i) {
      const auto i0 = static_cast<std::ptrdiff_t>(_chrom_bin_offsets[i - 1]);
      const auto i1 = static_cast<std::ptrdiff_t>(_chrom_bin_offsets[i]);
      const auto scale = std::sqrt(_scale[i - 1]);
      std::transform(_biases.begin() + i0, _biases.begin() + i1, biases.begin() + i0,
                     [&](const auto n) {
                       return n == 0 ? std::numeric_limits<double>::quiet_NaN() : n / scale;
                     });
    }
  }
  return biases;
}

inline void ICE::mask_cis_interactions(nonstd::span<const std::size_t> bin1_ids,
                                       nonstd::span<const std::size_t> bin2_ids,
                                       nonstd::span<double> counts,
                                       nonstd::span<const std::size_t> chrom_bin_offsets) {
  std::size_t j = 0;
  for (std::size_t i = 1; i < chrom_bin_offsets.size(); ++i) {
    const auto last_bin_id = chrom_bin_offsets[i];

    while (j < counts.size()) {
      if (bin1_ids[j] < last_bin_id) {
        if (bin2_ids[j] < last_bin_id) {
          counts[j] = 0;
        }
      } else {
        break;
      }
      ++j;
    }
  }
}

inline std::vector<double> ICE::scale() const noexcept { return _scale; }
inline std::vector<double> ICE::variance() const noexcept { return _variance; }

}  // namespace hictk::balancing
