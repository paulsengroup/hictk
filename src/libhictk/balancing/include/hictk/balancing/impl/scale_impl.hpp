// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "hictk/balancing/common.hpp"
#include "hictk/balancing/sparse_matrix.hpp"
#include "hictk/balancing/vc.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"

namespace hictk::balancing {

template <typename File>
inline SCALE::SCALE(const File& f, Type type, const Params& params) {
  internal::check_storage_mode(f);
  internal::check_bin_type(f);

  switch (type) {
    case Type::cis: {
      auto res = compute_cis(f, params);
      _chrom_offsets = std::move(res.offsets);
      _biases = std::move(res.weights);
      _scale = std::move(res.scales);
      return;
    }
    case Type::trans: {
      auto res = compute_trans(f, params);
      _chrom_offsets = std::move(res.offsets);
      _biases = std::move(res.weights);
      _scale = std::move(res.scales);
      return;
    }
    case Type::gw: {
      auto res = compute_gw(f, params);
      _chrom_offsets = std::move(res.offsets);
      _biases = std::move(res.weights);
      _scale = std::move(res.scales);
    }
  }
}

template <typename PixelIt>
inline SCALE::SCALE(PixelIt first, PixelIt last, const hictk::BinTable& bins, const Params& params)
    : _biases(VC{first, last, bins}.get_weights().to_vector(balancing::Weights::Type::DIVISIVE)),
      _convergence_stats(ConvergenceStats{false, false, 1000, 0, 10.0 * (1.0 + params.tol)}),
      _tpool(params.threads > 1 ? std::make_unique<BS::thread_pool>(params.threads) : nullptr) {
  if (bins.type() == BinTable::Type::variable) {
    throw std::runtime_error(
        "balancing interactions referring to a table with variable bin size is not supported");
  }
  if (first == last) {
    std::fill(_biases.begin(), _biases.end(), 1.0);
    _scale.push_back(1.0);
    _chrom_offsets = bins.num_bin_prefix_sum();
    return;
  }

  std::transform(_biases.begin(), _biases.end(), _biases.begin(),
                 [](const auto n) { return std::sqrt(n); });

  _max_tot_iters = params.max_iters * 3;

  const auto offset = bins.num_bin_prefix_sum().front();
  const auto matrix = mask_bins_and_init_buffers(first, last, offset, params.max_percentile,
                                                 params.tmpfile, params.chunk_size);

  std::visit([&](const auto& m) { balance(m, bins, params); }, matrix);
}

template <typename Matrix>
inline void SCALE::balance(const Matrix& m, const BinTable& bins, const Params& params) {
  VectorOfAtomicDecimals column(size(), 9);
  VectorOfAtomicDecimals row(size(), 9);

  m.multiply(row, _one, _tpool.get());
  row.multiply(_biases);

  auto dr = _biases;
  auto dc = _biases;
  auto current = _biases;

  std::vector<double> b_conv(size(), 0);
  std::vector<double> b0(size(), 0);
  std::vector<bool> bad_conv(size(), false);
  _ber_conv = 10.0;

  for (_iter = 0, _tot_iter = 0; _convergence_stats.error > params.tol &&
                                 _iter < params.max_iters && _tot_iter < _max_tot_iters;
       ++_iter, ++_tot_iter) {
    update_weights(column, _bad, row, _z_target_vector, dr, m, _tpool.get());
    column.multiply(dc);

    update_weights(row, _bad, column, _z_target_vector, dc, m, _tpool.get());
    row.multiply(dr);

    geometric_mean(dr, dc, _biases1);
    const auto res = compute_convergence_error(_biases1, current, _bad, params.tol);
    _convergence_stats.error = res.first;
    const auto num_bad = res.second;

    b0 = current;
    current = _biases1;

    _error_queue_iter.push(_convergence_stats.error);
    if (_error_queue_iter.size() == 7) {
      _error_queue_iter.pop();
    }

    const auto frac_bad = static_cast<double>(num_bad) / static_cast<double>(_nnz_rows);

    SPDLOG_INFO(FMT_STRING("Iteration {}: {}"), _tot_iter, _convergence_stats.error);

    if (_convergence_stats.error < params.tol) {
      SPDLOG_DEBUG(FMT_STRING("handle_convergence"));
      const auto status = handle_convergenece(m, dr, dc, row);
      if (status == ControlFlow::break_loop) {
        break;
      }
      assert(status == ControlFlow::continue_loop);
      reset_iter();
      continue;
    }

    if (_iter <= 4) {
      continue;
    }

    // check whether convergence rate is satisfactory
    const auto err1 = _error_queue_iter.front();
    const auto err2 = _error_queue_iter.back();
    if (err2 * (1.0 + params.delta) < err1 && (_iter < params.max_iters)) {
      continue;
    }

    // handle divergence
    SPDLOG_DEBUG(FMT_STRING("handle_divergence"));
    _convergence_stats.diverged = true;
    _convergence_stats.low_divergence = static_cast<std::uint32_t>(_low_cutoff);
    const auto status =
        handle_diverged(m, b0, dr, dc, row, frac_bad, params.frac_bad_cutoff, params.tol);
    if (status == ControlFlow::break_loop) {
      break;
    }
    if (status == ControlFlow::continue_loop) {
      continue;
    }
  }

  m.multiply(column, _biases1, _tpool.get());
  const auto row_sum_error = compute_final_error(column, _biases1, _z_target_vector, _bad);

  if (_convergence_stats.error > params.tol) {
    SPDLOG_DEBUG(FMT_STRING("error > tol: {} > {}"), _convergence_stats.error, params.tol);
  }
  if (row_sum_error > params.max_row_sum_error) {
    SPDLOG_DEBUG(FMT_STRING("row_sum_error > params.max_row_sum_error: {} > {}"), row_sum_error,
                 params.max_row_sum_error);
  }
  if (_low_cutoff > _upper_bound) {
    SPDLOG_DEBUG(FMT_STRING("low_cutoff > upper_bound: {} > {}"), _low_cutoff, _upper_bound);
  }
  // convergence not achieved, return vector of nans
  if (_convergence_stats.error > params.tol || row_sum_error > params.max_row_sum_error ||
      _low_cutoff > _upper_bound) {
    std::fill(_biases.begin(), _biases.end(), std::numeric_limits<double>::quiet_NaN());
    _scale.push_back(std::numeric_limits<double>::quiet_NaN());
    _chrom_offsets = bins.num_bin_prefix_sum();
    return;
  }

  // convergence achieved
  for (std::size_t i = 0; i < size(); ++i) {
    _biases[i] = _bad[i] ? std::numeric_limits<double>::quiet_NaN() : 1.0 / _biases1[i];
  }
  _scale.push_back(m.compute_scaling_factor_for_scale(_biases));
  _chrom_offsets = bins.num_bin_prefix_sum();
}

inline std::size_t SCALE::size() const noexcept { return _biases.size(); }

inline void SCALE::reset_iter() noexcept {
  _iter = 0;
  while (!_error_queue_iter.empty()) {
    _error_queue_iter.pop();
  }
}

inline balancing::Weights SCALE::get_weights(bool rescale) const {
  if (!rescale) {
    return {_biases, balancing::Weights::Type::DIVISIVE};
  }

  std::vector<double> biases(_biases.size());
  std::uint64_t chrom_id = 0;
  for (std::size_t i = 0; i < _biases.size(); ++i) {
    if (i >= _chrom_offsets[chrom_id + 1]) {
      chrom_id++;
    }
    biases[i] = _biases[i] * _scale[chrom_id];
  }

  return {biases, balancing::Weights::Type::DIVISIVE};
}

inline const std::vector<double>& SCALE::get_scale() const noexcept { return _scale; }

template <typename File>
inline auto SCALE::compute_cis(const File& f,
                               const hictk::balancing::SCALE::Params& params) -> Result {
  std::vector<std::uint64_t> offsets{};
  std::vector<double> scales{};
  std::vector<double> weights{};
  for (const Chromosome& chrom : f.chromosomes()) {
    if (chrom.is_all()) {
      continue;
    }
    const auto sel = f.fetch(chrom.name());
    const SCALE scale{sel.template begin<double>(), sel.template end<double>(),
                      f.bins().subset(chrom), params};

    offsets.push_back(f.bins().subset(chrom).num_bin_prefix_sum().front());

    const auto chrom_weights = scale.get_weights(false)(balancing::Weights::Type::DIVISIVE);
    scales.push_back(scale.get_scale().front());
    weights.insert(weights.end(), chrom_weights.begin(), chrom_weights.end());
  }

  offsets.push_back(f.bins().size());

  return {offsets, scales, weights};
}

template <typename File>
inline auto SCALE::compute_trans(const File& f, const Params& params) -> Result {
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
  const SCALE scale{sel.begin(), sel.end(), f.bins(), params};

  return {{0, f.bins().size()},
          scale.get_scale(),
          scale.get_weights(false).to_vector(balancing::Weights::Type::DIVISIVE)};
}

template <typename File>
inline auto SCALE::compute_gw(const File& f, const Params& params) -> Result {
  const auto sel = f.fetch();
  const SCALE scale{sel.template begin<double>(), sel.template end<double>(), f.bins(), params};

  return {{0, f.bins().size()},
          scale.get_scale(),
          scale.get_weights(false).to_vector(balancing::Weights::Type::DIVISIVE)};
}

template <typename Matrix>
inline void SCALE::update_weights(VectorOfAtomicDecimals& buffer, const std::vector<bool>& bad,
                                  VectorOfAtomicDecimals& weights,
                                  const std::vector<double>& target, std::vector<double>& d_vector,
                                  const Matrix& m, BS::thread_pool* tpool) noexcept {
  assert(buffer.size() == bad.size());
  assert(buffer.size() == weights.size());
  assert(buffer.size() == target.size());
  assert(buffer.size() == d_vector.size());

  for (std::size_t i = 0; i < weights.size(); ++i) {
    if (bad[i]) {
      weights.set(i, 1.0);
    }
    d_vector[i] *= target[i] / weights[i];
  }

  m.multiply(buffer, d_vector, tpool);
}

inline void SCALE::geometric_mean(const std::vector<double>& v1, const std::vector<double>& v2,
                                  std::vector<double>& vout) noexcept {
  assert(v1.size() == v2.size());
  assert(vout.size() == v1.size());

  for (std::size_t i = 0; i < v1.size(); ++i) {
    vout[i] = std::sqrt(v1[i] * v2[i]);
  }
}

inline std::pair<double, std::uint64_t> SCALE::compute_convergence_error(
    const std::vector<double>& biases, const std::vector<double>& current,
    const std::vector<bool>& bad, double tolerance) noexcept {
  assert(biases.size() == current.size());
  assert(biases.size() == bad.size());

  double error = 0;
  std::uint64_t num_fail = 0;
  for (std::size_t i = 0; i < biases.size(); ++i) {
    if (bad[i]) {
      continue;
    }
    const auto rel_err = std::abs((biases[i] - current[i]) / (biases[i] + current[i]));
    error = std::max(rel_err, error);
    num_fail += rel_err > tolerance;
  }

  return std::make_pair(error, num_fail);
}

inline double SCALE::compute_final_error(const VectorOfAtomicDecimals& col,
                                         const std::vector<double>& scale,
                                         const std::vector<double>& target,
                                         const std::vector<bool>& bad) noexcept {
  assert(col.size() == scale.size());
  assert(col.size() == target.size());
  assert(col.size() == bad.size());

  double error = 0.0;

  for (std::size_t i = 0; i < col.size(); ++i) {
    if (bad[i]) {
      continue;
    }
    const auto err1 = std::abs(col[i] * scale[i] - target[i]);
    error = std::max(error, err1);
  }

  return error;
}

inline void SCALE::multiply(std::vector<double>& v1, const std::vector<double>& v2) noexcept {
  assert(v1.size() == v2.size());
  for (std::size_t i = 0; i < v1.size(); ++i) {
    v1[i] *= v2[i];
  }
}

template <typename PixelIt>
inline std::variant<SparseMatrix, SparseMatrixChunked> SCALE::mask_bins_and_init_buffers(
    PixelIt first, PixelIt last, std::size_t offset, double max_percentile,
    const std::filesystem::path& tmpfile, std::size_t chunk_size) {
  assert(_bad.empty());
  assert(_one.empty());
  assert(_z_target_vector.empty());
  assert(_row_wise_nnz.empty());
  assert(_biases1.empty());

  // init buffers
  _bad.resize(size(), false);
  _one.resize(size(), 1);
  _z_target_vector.resize(size(), 1);
  _row_wise_nnz.resize(size(), 0);
  _biases1.resize(size(), 0);

  std::variant<SparseMatrix, SparseMatrixChunked> matrix{SparseMatrix{}};

  if (!tmpfile.empty()) {
    matrix = SparseMatrixChunked(tmpfile, chunk_size);
  }

  std::visit(
      [&](auto& m) {
        std::for_each(first, last, [&](const auto& p) {
          const auto bin1_id = p.bin1_id - offset;
          const auto bin2_id = p.bin2_id - offset;
          _row_wise_nnz[bin1_id]++;
          if (bin1_id != bin2_id) {
            _row_wise_nnz[bin2_id]++;
          }

          m.push_back(bin1_id, bin2_id, p.count);
        });
        m.finalize();
      },
      matrix);

  // compute the number of non-zero rows
  // we are sorting the vector of nnz because anyways we will need that in a later stage
  std::vector<std::uint64_t> row_wise_nnz_sorted{};
  std::copy_if(_row_wise_nnz.begin(), _row_wise_nnz.end(), std::back_inserter(row_wise_nnz_sorted),
               [&](const auto n) { return n != 0; });
  std::sort(row_wise_nnz_sorted.begin(), row_wise_nnz_sorted.end());
  _nnz_rows = row_wise_nnz_sorted.size() -
              static_cast<std::size_t>(std::distance(
                  row_wise_nnz_sorted.begin(),
                  std::upper_bound(row_wise_nnz_sorted.begin(), row_wise_nnz_sorted.end(), 0)));

  // compute the max number of nnz that can cause a row to be masked
  const auto upper_bound_idx =
      static_cast<std::size_t>(max_percentile * static_cast<double>(_nnz_rows) / 100.0);
  _upper_bound = row_wise_nnz_sorted[upper_bound_idx];

  // mask bad bins
  _low_cutoff = 1;
  for (std::size_t i = 0; i < _row_wise_nnz.size(); ++i) {
    if (_row_wise_nnz[i] < _low_cutoff) {
      _bad[i] = true;
      _one[i] = 0;
      _z_target_vector[i] = 0;
    }
  }
  return matrix;
}

template <typename Matrix>
inline auto SCALE::handle_convergenece(const Matrix& m, std::vector<double>& dr,
                                       std::vector<double>& dc,
                                       VectorOfAtomicDecimals& row) -> ControlFlow {
  _yes = true;
  if (_low_cutoff == 1) {
    SPDLOG_DEBUG(FMT_STRING("low cutoff"));
    return ControlFlow::break_loop;
  }
  SPDLOG_DEBUG(FMT_STRING("non low cutoff"));
  _convergence_stats.converged = true;
  _b_conv = _biases1;
  _bad_conv = _bad;
  _ber_conv = _convergence_stats.error;
  _convergence_stats.low_convergence = static_cast<std::uint32_t>(_low_cutoff);

  if (_convergence_stats.diverged) {
    if (_convergence_stats.low_convergence - _convergence_stats.low_divergence <= 1) {
      return ControlFlow::break_loop;
    }
    _low_cutoff = (_convergence_stats.low_convergence + _convergence_stats.low_divergence) / 2;
  } else {
    _low_cutoff = _convergence_stats.low_convergence / 2;
  }

  for (std::size_t i = 0; i < size(); ++i) {
    if (_row_wise_nnz[i] < _low_cutoff) {
      _bad[i] = true;
      _one[i] = 0;
    } else {
      _bad[i] = false;
      _one[i] = 1;
    }
  }

  _convergence_stats.error = 10.0;
  reset_iter();
  std::transform(_bad.begin(), _bad.end(), dr.begin(), [&](const auto b) { return !b; });
  dc = dr;

  m.multiply(row, dc, _tpool.get());
  row.multiply(dr);
  return ControlFlow::continue_loop;
}

template <typename Matrix>
inline auto SCALE::handle_almost_converged(const Matrix& m, const std::vector<double>& b0,
                                           std::vector<double>& dr, std::vector<double>& dc,
                                           VectorOfAtomicDecimals& row,
                                           double tolerance) -> ControlFlow {
  for (std::size_t i = 0; i < size(); ++i) {
    if (_bad[i]) {
      continue;
    }
    const auto rel_err = std::abs((_biases1[i] - b0[i]) / (_biases1[i] + b0[i]));
    if (rel_err > tolerance) {
      _bad[i] = true;
      _one[i] = 0;
    }
  }
  _yes = false;
  _convergence_stats.error = 10.0;
  reset_iter();
  std::transform(_bad.begin(), _bad.end(), dr.begin(), [&](const auto b) { return !b; });
  dc = dr;

  m.multiply(row, dc, _tpool.get());
  row.multiply(dr);

  if (_low_cutoff > _upper_bound) {
    return ControlFlow::break_loop;
  }
  if (_tot_iter > _max_tot_iters) {
    return ControlFlow::break_loop;
  }
  return ControlFlow::continue_loop;
}

template <typename Matrix>
inline auto SCALE::handle_diverged(const Matrix& m, const std::vector<double>& b0,
                                   std::vector<double>& dr, std::vector<double>& dc,
                                   VectorOfAtomicDecimals& row, double frac_bad,
                                   double frac_bad_cutoff, double tolerance) -> ControlFlow {
  const auto almost_converged = frac_bad < frac_bad_cutoff && _yes;
  if (_convergence_stats.converged) {
    if (_convergence_stats.low_convergence - _convergence_stats.low_divergence <= 1) {
      _biases1 = _b_conv;
      _bad = _bad_conv;
      _convergence_stats.error = _ber_conv;
      return ControlFlow::break_loop;
    }
    if (almost_converged) {
      const auto status = handle_almost_converged(m, b0, dr, dc, row, tolerance);
      if (status == ControlFlow::continue_loop) {
        return ControlFlow::continue_loop;
      }
      assert(status == ControlFlow::break_loop);
      return ControlFlow::break_loop;
    }
    _low_cutoff = (_convergence_stats.low_divergence + _convergence_stats.low_convergence) / 2;
    _yes = true;
  } else if (almost_converged) {
    const auto status = handle_almost_converged(m, b0, dr, dc, row, tolerance);
    if (status == ControlFlow::continue_loop) {
      return ControlFlow::continue_loop;
    }
    assert(status == ControlFlow::break_loop);
    return ControlFlow::break_loop;
  } else {
    _low_cutoff *= 2;
    _yes = true;
  }

  for (std::size_t i = 0; i < size(); ++i) {
    _bad[i] = _row_wise_nnz[i] < _low_cutoff;
    _one[i] = !_bad[i];
  }
  _convergence_stats.error = 10.0;
  reset_iter();

  dr = _one;
  dc = _one;
  m.multiply(row, dc, _tpool.get());
  row.multiply(dr);

  if (_low_cutoff > _upper_bound) {
    return ControlFlow::break_loop;
  }
  if (_tot_iter > _max_tot_iters) {
    return ControlFlow::break_loop;
  }
  return ControlFlow::continue_loop;
}

template <typename PixelIt>
inline std::variant<SparseMatrix, SparseMatrixChunked> SCALE::init_matrix(
    PixelIt first, PixelIt last, std::size_t offset, const std::filesystem::path& tmpfile,
    std::size_t chunk_size) {
  std::variant<SparseMatrix, SparseMatrixChunked> matrix{SparseMatrix{}};

  if (!tmpfile.empty()) {
    matrix = SparseMatrixChunked(tmpfile, chunk_size);
  }
  std::visit(
      [&](auto& m) {
        std::for_each(first, last, [&](const auto& p) {
          m.push_back(p.bin1_id - offset, p.bin2_id - offset, p.count);
        });
        m.finalize();
      },
      matrix);

  return matrix;
}

inline VC::Type SCALE::map_type_to_vc(Type type) noexcept {
  switch (type) {
    case Type::cis:
      return VC::Type::cis;
    case Type::trans:
      return VC::Type::trans;
    case Type::gw:
      return VC::Type::gw;
  }
}

}  // namespace hictk::balancing
