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
    : _chrom_offsets(f.bins().num_bin_prefix_sum()),
      _biases(f.bins().size(), 1.0),
      _variance(f.chromosomes().size(), 0),
      _scale(f.chromosomes().size(), std::numeric_limits<double>::quiet_NaN()) {
  std::unique_ptr<BS::thread_pool> tpool{};
  if (params.threads != 1) {
    tpool = std::make_unique<BS::thread_pool>(params.threads);
  }
  if (params.tmpfile.empty()) {
    balance_in_memory(f, type, params.tol, params.max_iters, params.num_masked_diags,
                      params.min_nnz, params.min_count, params.mad_max, tpool.get());
  } else {
    balance_chunked(f, type, params.tol, params.max_iters, params.num_masked_diags, params.min_nnz,
                    params.min_count, params.mad_max, params.tmpfile, params.chunk_size,
                    tpool.get());
  }
}

template <typename File>
inline void ICE::balance_in_memory(const File& f, Type type, double tol, std::size_t max_iters,
                                   std::size_t num_masked_diags, std::size_t min_nnz,
                                   std::size_t min_count, double mad_max, BS::thread_pool* tpool) {
  auto matrix = construct_sparse_matrix(f, type, num_masked_diags);

  initialize_biases(matrix, _biases, _chrom_offsets, min_nnz, min_count, mad_max, tpool);

  if (type == Type::gw) {
    return balance_gw(matrix, max_iters, tol, tpool);
  }

  if (type == Type::trans) {
    matrix.clear(true);
    matrix = construct_sparse_matrix_trans(f, num_masked_diags);
    return balance_trans(matrix, f.bins(), max_iters, tol, tpool);
  }

  assert(type == Type::cis);
  matrix.clear(true);
  for (std::uint32_t i = 0; i < f.chromosomes().size(); ++i) {
    const Chromosome& chrom = f.chromosomes().at(i);
    if (chrom.is_all()) {
      continue;
    }
    matrix = construct_sparse_matrix_cis(f, chrom, _chrom_offsets[i], num_masked_diags);
    balance_cis(matrix, chrom, max_iters, tol, tpool);
  }
}

template <typename File>
inline void ICE::balance_chunked(const File& f, Type type, double tol, std::size_t max_iters,
                                 std::size_t num_masked_diags, std::size_t min_nnz,
                                 std::size_t min_count, double mad_max,
                                 const std::filesystem::path& tmpfile, std::size_t chunk_size,
                                 BS::thread_pool* tpool) {
  auto matrix = construct_sparse_matrix_chunked(f, type, num_masked_diags, tmpfile, chunk_size);

  initialize_biases(matrix, _biases, _chrom_offsets, min_nnz, min_count, mad_max, tpool);

  if (type == Type::gw) {
    return balance_gw(matrix, max_iters, tol, tpool);
  }

  if (type == Type::trans) {
    matrix.clear(true);
    matrix = construct_sparse_matrix_chunked_trans(f, num_masked_diags, tmpfile, chunk_size);
    return balance_trans(matrix, f.bins(), max_iters, tol, tpool);
  }

  assert(type == Type::cis);
  matrix.clear(true);
  for (std::uint32_t i = 0; i < f.chromosomes().size(); ++i) {
    const Chromosome& chrom = f.chromosomes().at(i);
    if (chrom.is_all()) {
      continue;
    }
    matrix = construct_sparse_matrix_chunked_cis(f, chrom, _chrom_offsets[i], num_masked_diags,
                                                 tmpfile, chunk_size);
    balance_cis(matrix, chrom, max_iters, tol, tpool);
  }
}

template <typename MatrixT>
inline void ICE::balance_gw(const MatrixT& matrix, std::size_t max_iters, double tol,
                            BS::thread_pool* tpool) {
  _variance.resize(1, 0);
  _scale.resize(1, std::numeric_limits<double>::quiet_NaN());

  MargsVector marg(_biases.size());
  for (std::size_t i = 0; i < max_iters; ++i) {
    const auto res = inner_loop(matrix, _biases, marg, {}, tpool);
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
                               double tol, BS::thread_pool* tpool) {
  _variance.resize(1, 0);
  _scale.resize(1, std::numeric_limits<double>::quiet_NaN());
  const auto weights = compute_weights_from_chromosome_sizes(bins, _chrom_offsets);

  MargsVector marg(_biases.size());
  for (std::size_t i = 0; i < max_iters; ++i) {
    const auto res = inner_loop(matrix, _biases, marg, weights, tpool);
    SPDLOG_INFO(FMT_STRING("Iteration {}: {}"), i + 1, res.variance);
    _variance[0] = res.variance;
    _scale[0] = res.scale;
    if (res.variance < tol) {
      return;
    }
  }
}

template <typename MatrixT>
inline void ICE::balance_cis(const MatrixT& matrix, const Chromosome& chrom, std::size_t max_iters,
                             double tol, BS::thread_pool* tpool) {
  const auto i0 = _chrom_offsets[chrom.id()];
  const auto i1 = _chrom_offsets[chrom.id() + 1];
  auto biases_ = nonstd::span(_biases).subspan(i0, i1 - i0);

  MargsVector marg(biases_.size());
  for (std::size_t k = 0; k < max_iters; ++k) {
    const auto res = inner_loop(matrix, biases_, marg, {}, tpool);
    SPDLOG_INFO(FMT_STRING("[{}] iteration {}: {}"), chrom.name(), k + 1, res.variance);
    _variance[chrom.id()] = res.variance;
    _scale[chrom.id()] = res.scale;

    if (res.variance < tol) {
      break;
    }
  }
}

template <typename File>
auto ICE::construct_sparse_matrix(const File& f, Type type, std::size_t num_masked_diags)
    -> SparseMatrix {
  SPDLOG_INFO(FMT_STRING("Reading interactions into memory..."));
  if (type == Type::cis) {
    return construct_sparse_matrix_cis(f, num_masked_diags);
  }
  return construct_sparse_matrix_gw(f, num_masked_diags);
}

template <typename File>
inline auto ICE::construct_sparse_matrix_gw(const File& f, std::size_t num_masked_diags)
    -> SparseMatrix {
  SparseMatrix m{};

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
[[nodiscard]] inline auto ICE::construct_sparse_matrix_cis(const File& f, const Chromosome& chrom,
                                                           std::size_t bin_offset,
                                                           std::size_t num_masked_diags)
    -> SparseMatrix {
  SparseMatrix m{};

  const auto sel = f.fetch(chrom.name());
  std::for_each(sel.template begin<double>(), sel.template end<double>(),
                [&](const ThinPixel<double>& p) {
                  if (p.bin2_id - p.bin1_id >= num_masked_diags) {
                    m.push_back(p.bin1_id, p.bin2_id, p.count, bin_offset);
                  }
                });
  m.finalize();

  return m;
}

template <typename File>
[[nodiscard]] inline auto ICE::construct_sparse_matrix_cis(const File& f,
                                                           std::size_t num_masked_diags)
    -> SparseMatrix {
  SparseMatrix m{};

  for (const auto& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    auto sel = f.fetch(chrom.name());
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

  SparseMatrix m{};
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
  SPDLOG_INFO(FMT_STRING("Writing interactions to temporary file {}..."), tmpfile);
  if (type == Type::cis) {
    return construct_sparse_matrix_chunked_cis(f, num_masked_diags, tmpfile, chunk_size);
  }
  return construct_sparse_matrix_chunked_gw(f, num_masked_diags, tmpfile, chunk_size);
}

template <typename File>
inline auto ICE::construct_sparse_matrix_chunked_gw(const File& f, std::size_t num_masked_diags,
                                                    const std::filesystem::path& tmpfile,
                                                    std::size_t chunk_size) -> SparseMatrixChunked {
  SparseMatrixChunked m(tmpfile, chunk_size);

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
inline auto ICE::construct_sparse_matrix_chunked_cis(const File& f, const Chromosome& chrom,
                                                     std::size_t bin_offset,
                                                     std::size_t num_masked_diags,
                                                     const std::filesystem::path& tmpfile,
                                                     std::size_t chunk_size)
    -> SparseMatrixChunked {
  SparseMatrixChunked m(tmpfile, chunk_size);

  const auto sel = f.fetch(chrom.name());
  std::for_each(sel.template begin<double>(), sel.template end<double>(),
                [&](const ThinPixel<double>& p) {
                  if (p.bin2_id - p.bin1_id >= num_masked_diags) {
                    m.push_back(p.bin1_id, p.bin2_id, p.count, bin_offset);
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
  SparseMatrixChunked m(tmpfile, chunk_size);

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

  SparseMatrixChunked m(tmpfile, chunk_size);
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
inline void ICE::min_nnz_filtering(MargsVector& marg, const MatrixT& matrix,
                                   nonstd::span<double> biases, std::size_t min_nnz,
                                   BS::thread_pool* tpool) {
  matrix.marginalize_nnz(marg, tpool);
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (marg()[i] < static_cast<double>(min_nnz)) {
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
inline auto ICE::inner_loop(const MatrixT& matrix, nonstd::span<double> biases, MargsVector& marg,
                            nonstd::span<const double> weights, BS::thread_pool* tpool) -> Result {
  if (matrix.empty()) {
    std::fill(biases.begin(), biases.end(), std::numeric_limits<double>::quiet_NaN());
    return {std::numeric_limits<double>::quiet_NaN(), 0.0};
  }

  marg.resize(biases.size());
  matrix.times_outer_product_marg(marg, biases, weights, tpool);
  const auto [marg_sum, nnz_marg] = aggregate_marg(marg(), tpool);

  if (nnz_marg == 0) {
    std::fill(biases.begin(), biases.end(), std::numeric_limits<double>::quiet_NaN());
    return {std::numeric_limits<double>::quiet_NaN(), 0.0};
  }

  const auto avg_nzmarg = (marg_sum / static_cast<double>(nnz_marg));
  update_biases(marg(), biases, avg_nzmarg, tpool);

  const auto ssq_nzmarg = compute_ssq_nzmarg(marg(), avg_nzmarg, tpool);
  const auto var_nzmarg = ssq_nzmarg / static_cast<double>(nnz_marg - 1);

  return {avg_nzmarg, var_nzmarg};
}

inline std::pair<double, std::size_t> ICE::aggregate_marg(nonstd::span<const double> marg,
                                                          BS::thread_pool* tpool) {
  double marg_sum = 0.0;
  std::size_t nnz_marg{};

  std::mutex mtx{};

  auto aggregate_marg_impl = [&](std::size_t istart, std::size_t iend) {
    double marg_sum_ = 0.0;
    std::size_t nnz_marg_{};

    for (auto i = istart; i < iend; ++i) {
      marg_sum_ += marg[i];
      nnz_marg_ += marg[i] != 0;
    }

    [[maybe_unused]] const std::scoped_lock lck(mtx);
    marg_sum += marg_sum_;
    nnz_marg += nnz_marg_;
  };

  if (marg.size() < 10'000 || !tpool) {
    aggregate_marg_impl(0, marg.size());
    return std::make_pair(marg_sum, nnz_marg);
  }

  tpool->push_loop(0, marg.size(), aggregate_marg_impl);
  tpool->wait_for_tasks();

  return std::make_pair(marg_sum, nnz_marg);
}

inline void ICE::update_biases(nonstd::span<const double> marg, nonstd::span<double> biases,
                               double avg_nzmarg, BS::thread_pool* tpool) {
  auto update_biases_impl = [&](std::size_t istart, std::size_t iend) {
    for (auto i = istart; i < iend; ++i) {
      const auto n = marg[i] / avg_nzmarg;
      if (n != 0) {
        biases[i] /= n;
      }
    }
  };

  if (marg.size() < 10'000 || !tpool) {
    return update_biases_impl(0, marg.size());
  }

  tpool->push_loop(0, marg.size(), update_biases_impl);
  tpool->wait_for_tasks();
}

inline double ICE::compute_ssq_nzmarg(nonstd::span<const double> marg, double avg_nzmarg,
                                      BS::thread_pool* tpool) {
  std::mutex mtx{};
  double ssq_nzmarg = 0;

  auto compute_ssq_nzmarg_impl = [&](std::size_t istart, std::size_t iend) {
    double ssq_nzmarg_ = 0.0;
    for (auto i = istart; i < iend; ++i) {
      const auto& n = marg[i];
      if (n != 0) {
        ssq_nzmarg_ += std::pow(n - avg_nzmarg, 2);
      }
    }
    [[maybe_unused]] const std::scoped_lock lck(mtx);
    ssq_nzmarg += ssq_nzmarg_;
  };

  if (marg.size() < 10'000 || !tpool) {
    compute_ssq_nzmarg_impl(0, marg.size());
    return ssq_nzmarg;
  }

  tpool->push_loop(0, marg.size(), compute_ssq_nzmarg_impl);
  tpool->wait_for_tasks();
  return ssq_nzmarg;
}

template <typename MatrixT>
inline void ICE::initialize_biases(const MatrixT& matrix, nonstd::span<double> biases,
                                   nonstd::span<const std::size_t> chrom_bin_offsets,
                                   std::size_t min_nnz, std::size_t min_count, double mad_max,
                                   BS::thread_pool* tpool) {
  if (min_nnz == 0 && min_count == 0 && mad_max == 0) {
    return;
  }

  SPDLOG_INFO(FMT_STRING("Initializing bias vector..."));
  MargsVector marg(biases.size());
  if (min_nnz != 0) {
    SPDLOG_INFO(FMT_STRING("Masking rows with fewer than {} nnz entries..."), min_nnz);
    min_nnz_filtering(marg, matrix, biases, min_nnz, tpool);
  }

  if (min_count != 0 || mad_max != 0) {
    matrix.marginalize(marg, tpool);
  }

  if (min_count != 0) {
    SPDLOG_INFO(FMT_STRING("Masking rows with fewer than {} interactions..."), min_count);
    min_count_filtering(biases, min_count, marg());
  }

  if (mad_max != 0) {
    SPDLOG_INFO(FMT_STRING("Masking rows using mad_max={}..."), mad_max);
    mad_max_filtering(chrom_bin_offsets, biases, marg(), mad_max);
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
