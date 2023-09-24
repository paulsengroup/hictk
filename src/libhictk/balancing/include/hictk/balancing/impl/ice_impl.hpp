// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>
#include <zstd.h>

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

template <typename File>
inline ICE::ICE(const File& f, Type type, const Params& params)
    : _chrom_offsets(f.bins().num_bin_prefix_sum()), _biases(f.bins().size(), 1.0) {
  if (params.tmpfile.empty()) {
    balance_in_memory(f, type, params.tol, params.max_iters, params.num_masked_diags,
                      params.min_nnz, params.min_count, params.mad_max);
  } else {
    balance_chunked(f, type, params.tol, params.max_iters, params.num_masked_diags, params.min_nnz,
                    params.min_count, params.mad_max, params.tmpfile, params.chunk_size);
  }
}

template <typename File>
inline void ICE::balance_in_memory(const File& f, Type type, double tol, std::size_t max_iters,
                                   std::size_t num_masked_diags, std::size_t min_nnz,
                                   std::size_t min_count, double mad_max) {
  auto matrix = construct_sparse_matrix(f, type, num_masked_diags);

  initialize_biases(matrix.view(), _biases, _chrom_offsets, min_nnz, min_count, mad_max);

  switch (type) {
    case Type::gw:
      balance_gw(matrix.view(), max_iters, tol);
      break;
    case Type::cis:
      balance_cis(matrix, f.bins(), max_iters, tol);
      break;
    case Type::trans:
      matrix = construct_sparse_matrix_trans(f, num_masked_diags);
      balance_trans(matrix.view(), f.bins(), max_iters, tol);
  }
}

template <typename File>
inline void ICE::balance_chunked(const File& f, Type type, double tol, std::size_t max_iters,
                                 std::size_t num_masked_diags, std::size_t min_nnz,
                                 std::size_t min_count, double mad_max,
                                 const std::filesystem::path& tmpfile, std::size_t chunk_size) {
  auto matrix = construct_sparse_matrix_chunked(f, type, num_masked_diags, tmpfile, chunk_size);

  initialize_biases(matrix.view(), _biases, _chrom_offsets, min_nnz, min_count, mad_max);

  switch (type) {
    case Type::gw:
      balance_gw(matrix.view(), max_iters, tol);
      break;
    case Type::cis:
      balance_cis(matrix, f.bins(), max_iters, tol);
      break;
    case Type::trans:
      balance_trans(
          construct_sparse_matrix_chunked_trans(f, num_masked_diags, tmpfile, chunk_size).view(),
          f.bins(), max_iters, tol);
  }
}

template <typename MatrixT>
inline void ICE::balance_gw(const MatrixT& matrix, std::size_t max_iters, double tol) {
  _variance.resize(1, 0);
  _scale.resize(1, std::numeric_limits<double>::quiet_NaN());

  for (std::size_t i = 0; i < max_iters; ++i) {
    const auto res = inner_loop(matrix, _biases);
    SPDLOG_INFO(FMT_STRING("Iteration {}: {}"), i + 1, res.variance);
    _variance[0] = res.variance;
    _scale[0] = res.scale;
    if (res.variance < tol) {
      return;
    }
  }
}

template <typename MatrixT>
inline void ICE::balance_trans(const MatrixT& matrix, const BinTable& bins, std::size_t max_iters,
                               double tol) {
  _variance.resize(1, 0);
  _scale.resize(1, std::numeric_limits<double>::quiet_NaN());
  const auto weights = compute_weights_from_chromosome_sizes(bins, _chrom_offsets);

  for (std::size_t i = 0; i < max_iters; ++i) {
    const auto res = inner_loop(matrix, _biases, weights);
    SPDLOG_INFO(FMT_STRING("Iteration {}: {}"), i + 1, res.variance);
    _variance[0] = res.variance;
    _scale[0] = res.scale;
    if (res.variance < tol) {
      return;
    }
  }
}

template <typename MatrixT>
inline void ICE::balance_cis(const MatrixT& matrix, const BinTable& bins, std::size_t max_iters,
                             double tol) {
  _variance.resize(bins.chromosomes().size(), 0);
  _scale.resize(bins.chromosomes().size(), std::numeric_limits<double>::quiet_NaN());

  std::vector<double> margs(_biases.size());

  for (const auto& chrom : bins.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto cis_matrix = matrix.subset(chrom.id());

    const auto j0 = _chrom_offsets[chrom.id()];
    const auto j1 = _chrom_offsets[chrom.id() + 1];

    auto biases_ = nonstd::span(_biases).subspan(j0, j1 - j0);

    for (std::size_t k = 0; k < max_iters; ++k) {
      const auto res = inner_loop(cis_matrix, biases_);
      SPDLOG_INFO(FMT_STRING("[{}] iteration {}: {}"), chrom.name(), k + 1, res.variance);
      _variance[chrom.id()] = res.variance;
      _scale[chrom.id()] = res.scale;

      if (res.variance < tol) {
        break;
      }
    }
  }
}

template <typename File>
auto ICE::construct_sparse_matrix(const File& f, Type type, std::size_t num_masked_diags)
    -> SparseMatrix {
  if (type == Type::cis) {
    return construct_sparse_matrix_cis(f, num_masked_diags);
  }
  return construct_sparse_matrix_gw(f, num_masked_diags);
}

template <typename File>
inline auto ICE::construct_sparse_matrix_gw(const File& f, std::size_t num_masked_diags)
    -> SparseMatrix {
  SparseMatrix m(f.bins());

  const auto sel = f.fetch();
  std::for_each(sel.template begin<double>(), sel.template end<double>(), [&](const auto& p) {
    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
  });

  m.finalize();

  return m;
}

template <typename File>
[[nodiscard]] inline auto ICE::construct_sparse_matrix_cis(const File& f,
                                                           std::size_t num_masked_diags)
    -> SparseMatrix {
  SparseMatrix m(f.bins());

  for (const Chromosome& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto sel = f.fetch(chrom.name());
    std::for_each(sel.template begin<double>(), sel.template end<double>(),
                  [&](const ThinPixel<double>& p) {
                    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
                      m.push_back(p.bin1_id, p.bin2_id, p.count);
                    }
                  });
  }
  m.finalize();

  return m;
}

template <typename File>
[[nodiscard]] inline auto ICE::construct_sparse_matrix_trans(const File& f,
                                                             std::size_t num_masked_diags)
    -> SparseMatrix {
  using SelectorT = decltype(f.fetch("chr1", "chr2"));
  using PixelIt = decltype(f.fetch("chr1", "chr2").template begin<double>());

  std::vector<SelectorT> selectors{};
  for (const Chromosome& chrom1 : f.chromosomes()) {
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1.id() + 1; chrom2_id < f.chromosomes().size();
         ++chrom2_id) {
      const auto& chrom2 = f.chromosomes().at(chrom2_id);
      if (chrom2.is_all()) {
        continue;
      }

      selectors.emplace_back(f.fetch(chrom1.name(), chrom2.name()));
    }
    std::vector<PixelIt> heads{};
    std::vector<PixelIt> tails{};
    for (const auto& sel : selectors) {
      heads.emplace_back(sel.template begin<double>());
      tails.emplace_back(sel.template end<double>());
    }
  }

  std::vector<PixelIt> heads{};
  std::vector<PixelIt> tails{};
  for (const auto& sel : selectors) {
    heads.emplace_back(sel.template begin<double>());
    tails.emplace_back(sel.template end<double>());
  }

  internal::PixelMerger<PixelIt> merger{heads, tails};

  SparseMatrix m(f.bins());
  std::for_each(merger.begin(), merger.end(), [&](const ThinPixel<double>& p) {
    // TODO: this filtering step is wrong when done on trans matrices, as it will
    // remove the first and last few pixels from trans matrices of adjacent chromosomes.
    // Remove the filtering once this bug has been fixed in cooler
    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
  });

  m.shrink_to_fit();

  return m;
}

template <typename File>
auto ICE::construct_sparse_matrix_chunked(const File& f, Type type, std::size_t num_masked_diags,
                                          const std::filesystem::path& tmpfile,
                                          std::size_t chunk_size) -> SparseMatrixChunked {
  if (type == Type::cis) {
    return construct_sparse_matrix_chunked_cis(f, num_masked_diags, tmpfile, chunk_size);
  }
  return construct_sparse_matrix_chunked_gw(f, num_masked_diags, tmpfile, chunk_size);
}

template <typename File>
inline auto ICE::construct_sparse_matrix_chunked_gw(const File& f, std::size_t num_masked_diags,
                                                    const std::filesystem::path& tmpfile,
                                                    std::size_t chunk_size) -> SparseMatrixChunked {
  SparseMatrixChunked m(f.bins(), tmpfile, chunk_size);

  const auto sel = f.fetch();
  std::for_each(sel.template begin<double>(), sel.template end<double>(), [&](const auto& p) {
    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
  });

  m.finalize();
  return m;
}

template <typename File>
inline auto ICE::construct_sparse_matrix_chunked_cis(const File& f, std::size_t num_masked_diags,
                                                     const std::filesystem::path& tmpfile,
                                                     std::size_t chunk_size)
    -> SparseMatrixChunked {
  SparseMatrixChunked m(f.bins(), tmpfile, chunk_size);

  for (const Chromosome& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto sel = f.fetch(chrom.name());
    std::for_each(sel.template begin<double>(), sel.template end<double>(),
                  [&](const ThinPixel<double>& p) {
                    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
                      m.push_back(p.bin1_id, p.bin2_id, p.count);
                    }
                  });
  }
  m.finalize();
  return m;
}

template <typename File>
inline auto ICE::construct_sparse_matrix_chunked_trans(const File& f, std::size_t num_masked_diags,
                                                       const std::filesystem::path& tmpfile,
                                                       std::size_t chunk_size)
    -> SparseMatrixChunked {
  using SelectorT = decltype(f.fetch("chr1", "chr2"));
  using PixelIt = decltype(f.fetch("chr1", "chr2").template begin<double>());

  std::vector<SelectorT> selectors{};
  for (const Chromosome& chrom1 : f.chromosomes()) {
    if (chrom1.is_all()) {
      continue;
    }
    for (std::uint32_t chrom2_id = chrom1.id() + 1; chrom2_id < f.chromosomes().size();
         ++chrom2_id) {
      const auto& chrom2 = f.chromosomes().at(chrom2_id);
      if (chrom2.is_all()) {
        continue;
      }

      selectors.emplace_back(f.fetch(chrom1.name(), chrom2.name()));
    }
  }

  std::vector<PixelIt> heads{};
  std::vector<PixelIt> tails{};
  for (const auto& sel : selectors) {
    heads.emplace_back(sel.template begin<double>());
    tails.emplace_back(sel.template end<double>());
  }

  internal::PixelMerger<PixelIt> merger{heads, tails};

  SparseMatrixChunked m(f.bins(), tmpfile, chunk_size);
  std::for_each(merger.begin(), merger.end(), [&](const ThinPixel<double>& p) {
    // TODO: this filtering step is wrong when done on trans matrices, as it will
    // remove the first and last few pixels from trans matrices of adjacent chromosomes.
    // Remove the filtering once this bug has been fixed in cooler
    if (p.bin2_id - p.bin1_id >= num_masked_diags) {
      m.push_back(p.bin1_id, p.bin2_id, p.count);
    }
  });

  m.finalize();
  return m;
}

template <typename MatrixT>
inline void ICE::min_nnz_filtering(const MatrixT& matrix, nonstd::span<double> biases,
                                   std::size_t min_nnz) {
  const auto& marg = matrix.marginalize_nnz();
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (marg[i] < static_cast<double>(min_nnz)) {
      biases[i] = 0;
    }
  }
}

inline void ICE::min_count_filtering(nonstd::span<double> biases, std::size_t min_count,
                                     nonstd::span<const double> marg) {
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (marg[i] < static_cast<double>(min_count)) {
      biases[i] = 0;
    }
  }
}

inline void ICE::mad_max_filtering(nonstd::span<const std::size_t> chrom_offsets,
                                   nonstd::span<double> biases, nonstd::span<double> marg,
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

template <typename MatrixT>
inline auto ICE::inner_loop(const MatrixT& matrix, nonstd::span<double> biases,
                            nonstd::span<const double> weights) -> Result {
  if (matrix.empty()) {
    std::fill(biases.begin(), biases.end(), std::numeric_limits<double>::quiet_NaN());
    return {std::numeric_limits<double>::quiet_NaN(), 0.0};
  }
  const auto& marg = matrix.times_outer_product_marg(biases, weights);

  double marg_sum = 0.0;
  std::size_t nnz_marg{};
  for (const auto& n : marg) {
    marg_sum += n;
    nnz_marg += n != 0;
  }

  if (nnz_marg == 0) {
    std::fill(biases.begin(), biases.end(), std::numeric_limits<double>::quiet_NaN());
    return {std::numeric_limits<double>::quiet_NaN(), 0.0};
  }

  const auto avg_nzmarg = (marg_sum / static_cast<double>(nnz_marg));
  for (std::size_t i = 0; i < biases.size(); ++i) {
    const auto n = marg[i] / avg_nzmarg;
    if (n != 0) {
      biases[i] /= n;
    }
  }

  double ssq_nzmarg = 0.0;
  for (const auto n : marg) {
    if (n != 0) {
      ssq_nzmarg += std::pow(n - avg_nzmarg, 2);
    }
  }
  const auto var_nzmarg = ssq_nzmarg / static_cast<double>(nnz_marg - 1);

  return {avg_nzmarg, var_nzmarg};
}

template <typename MatrixT>
inline void ICE::initialize_biases(const MatrixT& matrix, nonstd::span<double> biases,
                                   nonstd::span<const std::size_t> chrom_bin_offsets,
                                   std::size_t min_nnz, std::size_t min_count, double mad_max) {
  if (min_nnz != 0) {
    min_nnz_filtering(matrix, biases, min_nnz);
  }

  if (min_count != 0 || mad_max != 0) {
    matrix.marginalize();
  }
  if (min_count != 0) {
    min_count_filtering(biases, min_count, matrix.margs());
  }

  if (mad_max != 0) {
    auto margs = std::vector<double>{matrix.margs()};
    mad_max_filtering(chrom_bin_offsets, biases, margs, mad_max);
  }
}

inline std::vector<double> ICE::compute_weights_from_chromosome_sizes(
    const BinTable& bins, nonstd::span<std::size_t> chrom_bin_offsets) {
  std::vector<double> weights(bins.size());
  for (std::uint32_t i = 1; i < chrom_bin_offsets.size(); ++i) {
    const auto& chrom = bins.chromosomes().at(i - 1);
    if (chrom.is_all()) {
      continue;
    }
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

inline std::vector<double> ICE::get_weights(bool rescale) const {
  if (!rescale) {
    return _biases;
  }

  std::vector<double> biases(_biases.size());
  if (_scale.size() == 1) {
    const auto scale = std::sqrt(_scale[0]);
    std::transform(_biases.begin(), _biases.end(), biases.begin(), [&](const auto n) {
      return n == 0 ? std::numeric_limits<double>::quiet_NaN() : n / scale;
    });
  } else {
    for (std::size_t i = 1; i < _chrom_offsets.size(); ++i) {
      const auto i0 = static_cast<std::ptrdiff_t>(_chrom_offsets[i - 1]);
      const auto i1 = static_cast<std::ptrdiff_t>(_chrom_offsets[i]);
      const auto scale = std::sqrt(_scale[i - 1]);
      std::transform(_biases.begin() + i0, _biases.begin() + i1, biases.begin() + i0,
                     [&](const auto n) {
                       return n == 0 ? std::numeric_limits<double>::quiet_NaN() : n / scale;
                     });
    }
  }
  return biases;
}

inline std::vector<double> ICE::scale() const noexcept { return _scale; }
inline std::vector<double> ICE::variance() const noexcept { return _variance; }

}  // namespace hictk::balancing
