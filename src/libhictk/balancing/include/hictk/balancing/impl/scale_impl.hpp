// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <limits>
#include <utility>
#include <vector>

#include "hictk/balancing/vc.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"

namespace hictk::balancing {

template <typename File>
inline SCALE::SCALE(const File& f, Type type, const Params& params) {
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
    : _biases(VC{first, last, bins}.get_weights()) {
  if (first == last) {
    std::fill(_biases.begin(), _biases.end(), 1.0);
    _scale.push_back(1.0);
    _chrom_offsets = bins.num_bin_prefix_sum();
    return;
  }

  const auto offset = bins.num_bin_prefix_sum().front();
  const auto size = _biases.size();

  std::vector<bool> bad(size, false);

  std::vector<double> z_target_vector(size, 1);
  std::vector<double> calculated_vector_b(size, 0);
  std::vector<double> one(size, 1);

  const auto total_iters = params.max_iters * 3;

  std::vector<double> report_error_for_iteration{};
  std::vector<std::uint32_t> num_iters_for_all_iterations{};

  std::vector<std::size_t> row_wise_nnz(size, 0);
  std::for_each(first, last, [&](const auto& p) {
    row_wise_nnz[p.bin1_id - offset]++;
    if (p.bin1_id != p.bin2_id) {
      row_wise_nnz[p.bin2_id - offset]++;
    }
  });

  std::vector<std::size_t> row_wise_nnz_sorted{};
  std::copy_if(row_wise_nnz.begin(), row_wise_nnz.end(), std::back_inserter(row_wise_nnz_sorted),
               [&](const auto n) { return n != 0; });
  std::sort(row_wise_nnz_sorted.begin(), row_wise_nnz_sorted.end());
  const auto nnz_rows =
      row_wise_nnz_sorted.size() -
      static_cast<std::size_t>(std::distance(
          row_wise_nnz_sorted.begin(),
          std::upper_bound(row_wise_nnz_sorted.begin(), row_wise_nnz_sorted.end(), 0)));

  const auto upper_bound_idx =
      static_cast<std::size_t>(params.max_percentile * static_cast<double>(nnz_rows) / 100.0);
  const auto upper_bound = row_wise_nnz_sorted[upper_bound_idx];

  std::size_t low_cutoff = 1;

  for (std::size_t i = 0; i < row_wise_nnz.size(); ++i) {
    if (row_wise_nnz[i] < low_cutoff) {
      bad[i] = true;
      one[i] = 0;
      z_target_vector[i] = 0;
    }
  }

  // TODO avoid reading all interactions into memory
  const std::vector<ThinPixel<double>> pixels(first, last);

  auto row = matrix_vect_mult(pixels, one, offset, 0);
  multiply(row, _biases);

  for (std::size_t i = 0; i < size; ++i) {
    one[i] = 1 - bad[i];
  }

  auto dr = _biases;
  auto dc = _biases;
  auto current = _biases;

  bool conv = false;
  bool div = false;
  std::uint32_t low_conv = 1000;
  std::uint32_t low_div = 0;

  std::vector<double> b_conv(size, 0);
  std::vector<double> b0(size, 0);
  std::vector<bool> bad_conv(size, false);
  auto ber_conv = 10.0;
  bool yes = true;

  auto converge_error = 10.0 * (1.0 + params.tol);

  for (std::size_t iter = 0, all_iters_i = 0, real_iters = 0;
       converge_error > params.tol && iter < params.max_iters && all_iters_i < total_iters;
       ++iter, ++all_iters_i, ++real_iters) {
    auto col = update_weights(bad, row, z_target_vector, dr, pixels, offset);

    for (std::size_t i = 0; i < col.size(); ++i) {
      col[i] *= dc[i];
    }

    row = update_weights(bad, col, z_target_vector, dc, pixels, offset);

    for (std::size_t i = 0; i < col.size(); ++i) {
      row[i] *= dr[i];
    }

    geometric_mean(dr, dc, calculated_vector_b);
    const auto res = compute_convergence_error(calculated_vector_b, current, bad, params.tol);
    converge_error = res.first;
    const auto num_bad = res.second;

    b0 = current;

    report_error_for_iteration.push_back(converge_error);
    num_iters_for_all_iterations.push_back(static_cast<std::uint32_t>(iter));

    current = calculated_vector_b;
    const auto frac_bad = static_cast<double>(num_bad) / static_cast<double>(nnz_rows);

    if (converge_error < params.tol) {
      yes = true;
      if (low_cutoff == 1) {
        break;
      }
      conv = true;
      b_conv = calculated_vector_b;
      bad_conv = bad;
      ber_conv = converge_error;
      low_conv = static_cast<std::uint32_t>(low_cutoff);

      if (div) {
        if (low_conv - low_div <= 1) {
          break;
        }
        low_cutoff = (low_conv + low_div) / 2;
      } else {
        low_cutoff = low_conv / 2;
      }

      for (std::size_t i = 0; i < size; ++i) {
        if (row_wise_nnz[i] < low_cutoff) {
          bad[i] = true;
          one[i] = 0;
        } else {
          bad[i] = false;
          one[i] = 1;
        }
      }

      converge_error = 10.0;
      iter = 0;
      std::transform(bad.begin(), bad.end(), dr.begin(), [&](const auto b) { return !b; });
      dc = dr;

      row = matrix_vect_mult(pixels, dc, offset);
      multiply(row, dr);
      continue;
    }

    if (iter <= 5) {
      continue;
    }

    const auto err2 = report_error_for_iteration[report_error_for_iteration.size() - 1];
    const auto err1 = report_error_for_iteration[report_error_for_iteration.size() - 6];

    if (err2 * (1.0 + params.delta) < err1 && (iter < params.max_iters)) {
      continue;
    }

    div = true;
    low_div = static_cast<std::uint32_t>(low_cutoff);
    if (conv) {
      if (low_conv - low_div <= 1) {
        calculated_vector_b = b_conv;
        bad = bad_conv;
        converge_error = ber_conv;
        break;

      } else if (frac_bad < params.erez && yes) {
        for (std::size_t i = 0; i < size; ++i) {
          if (bad[i]) {
            continue;
          }
          const auto rel_err =
              std::abs((calculated_vector_b[i] - b0[i]) / (calculated_vector_b[i] + b0[i]));
          if (rel_err > params.tol) {
            bad[i] = true;
            one[i] = 0;
          }
        }
        yes = false;
        converge_error = 10.0;
        iter = 0;
        std::transform(bad.begin(), bad.end(), dr.begin(), [&](const auto b) { return !b; });
        dc = dr;

        row = matrix_vect_mult(pixels, dc, offset);
        multiply(row, dr);

        if (low_cutoff > upper_bound) {
          break;
        }
        if (all_iters_i > total_iters) {
          break;
        }
        continue;
      } else {
        low_cutoff = (low_div + low_conv) / 2;
        yes = true;
      }
    } else if (frac_bad < params.erez && yes) {
      for (std::size_t i = 0; i < size; ++i) {
        if (bad[i]) {
          continue;
        }
        const auto rel_err =
            std::abs((calculated_vector_b[i] - b0[i]) / (calculated_vector_b[i] + b0[i]));
        if (rel_err > params.tol) {
          bad[i] = true;
          one[i] = 0;
        }
      }

      yes = false;
      converge_error = 10.0;
      iter = 0;

      std::transform(bad.begin(), bad.end(), dr.begin(), [&](const auto b) { return !b; });
      dc = dr;
      row = matrix_vect_mult(pixels, dc, offset);
      multiply(row, dr);

      if (low_cutoff > upper_bound) {
        break;
      }
      if (all_iters_i > total_iters) {
        break;
      }
      continue;
    } else {
      low_cutoff *= 2;
      yes = true;
    }

    for (std::size_t i = 0; i < size; ++i) {
      if (row_wise_nnz[i] < low_cutoff) {
        bad[i] = true;
        one[i] = 0;
      } else {
        bad[i] = false;
        one[i] = 1;
      }
    }
    converge_error = 10.0;
    iter = 0;

    std::transform(bad.begin(), bad.end(), dr.begin(), [&](const auto b) { return !b; });
    dc = dr;
    row = matrix_vect_mult(pixels, dc, offset);
    multiply(row, dr);

    if (low_cutoff > upper_bound) {
      break;
    }
    if (all_iters_i > total_iters) {
      break;
    }
  }

  const auto col = matrix_vect_mult(pixels, calculated_vector_b, offset);
  const auto row_sum_error = compute_final_error(col, calculated_vector_b, z_target_vector, bad);

  if (converge_error > params.tol || row_sum_error > params.max_row_sum_error ||
      low_cutoff > upper_bound) {
    std::fill(_biases.begin(), _biases.end(), std::numeric_limits<double>::quiet_NaN());
    _scale.push_back(std::numeric_limits<double>::quiet_NaN());
    _chrom_offsets = bins.num_bin_prefix_sum();
    return;
  }

  for (std::size_t i = 0; i < size; ++i) {
    if (bad[i]) {
      _biases[i] = std::numeric_limits<double>::quiet_NaN();
    } else {
      _biases[i] = 1.0 / calculated_vector_b[i];
    }
  }
  _scale.push_back(compute_scale(pixels, _biases, offset));
  _chrom_offsets = bins.num_bin_prefix_sum();
}

inline std::vector<double> SCALE::get_weights(bool rescale) const {
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

  return biases;
}

inline const std::vector<double>& SCALE::get_scale() const noexcept { return _scale; }

template <typename File>
inline auto SCALE::compute_cis(const File& f, const hictk::balancing::SCALE::Params& params)
    -> Result {
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

    const auto chrom_weights = scale.get_weights(false);
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

  return {{0, f.bins().size()}, scale.get_scale(), scale.get_weights(false)};
}

template <typename File>
inline auto SCALE::compute_gw(const File& f, const Params& params) -> Result {
  const auto sel = f.fetch();
  const SCALE scale{sel.template begin<double>(), sel.template end<double>(), f.bins(), params};

  return {{0, f.bins().size()}, scale.get_scale(), scale.get_weights(false)};
}
inline std::vector<double> SCALE::matrix_vect_mult(const std::vector<ThinPixel<double>>& pixels,
                                                   const std::vector<double>& cfx,
                                                   std::size_t offset, std::size_t i0,
                                                   std::size_t i1) {
  std::vector<double> v(cfx.size());
  matrix_vect_mult(pixels, cfx, v, offset, i0, i1);
  return v;
}

inline void SCALE::matrix_vect_mult(const std::vector<ThinPixel<double>>& pixels,
                                    const std::vector<double>& cfx, std::vector<double>& sum_vect,
                                    std::size_t offset, std::size_t i0, std::size_t i1) noexcept {
  std::fill(sum_vect.begin(), sum_vect.end(), 0);

  if (pixels.empty()) {
    return;
  }

  if (i1 == 0) {
    i1 = pixels.size();
  }

  assert(i1 <= pixels.size());
  assert(i0 < i1);

  for (std::size_t i = i0; i < i1; ++i) {
    const auto& p = pixels[i];

    assert(p.bin1_id - offset < sum_vect.size());
    assert(p.bin2_id - offset < sum_vect.size());

    const auto f = p.bin1_id == p.bin2_id ? 0.5 : 1.0;
    sum_vect[p.bin1_id - offset] += p.count * f * cfx[p.bin2_id - offset];
    sum_vect[p.bin2_id - offset] += p.count * f * cfx[p.bin1_id - offset];
  }
}

inline std::vector<double> SCALE::update_weights(const std::vector<bool>& bad,
                                                 std::vector<double>& weights,
                                                 const std::vector<double>& target,
                                                 std::vector<double>& d_vector,
                                                 const std::vector<ThinPixel<double>>& pixels,
                                                 std::size_t offset) noexcept {
  assert(bad.size() == weights.size());
  for (std::size_t i = 0; i < weights.size(); ++i) {
    if (bad[i] == true) {
      weights[i] = 1.0;
    }
  }

  for (std::size_t i = 0; i < weights.size(); ++i) {
    d_vector[i] *= target[i] / weights[i];
  }

  return matrix_vect_mult(pixels, d_vector, offset, 0);
}

inline void SCALE::geometric_mean(const std::vector<double>& v1, const std::vector<double>& v2,
                                  std::vector<double>& vout) noexcept {
  assert(v1.size() == v2.size());
  assert(vout.size() == v1.size());

  std::fill(vout.begin(), vout.end(), 0);

  for (std::size_t i = 0; i < v1.size(); ++i) {
    vout[i] = std::sqrt(v1[i] * v2[i]);
  }
}

inline std::pair<double, std::uint64_t> SCALE::compute_convergence_error(
    const std::vector<double>& calculated_vector_b, const std::vector<double>& current,
    const std::vector<bool>& bad, double tolerance) noexcept {
  assert(calculated_vector_b.size() == current.size());
  assert(calculated_vector_b.size() == bad.size());

  double error = 0;
  std::uint64_t num_fail = 0;
  for (std::size_t i = 0; i < calculated_vector_b.size(); ++i) {
    if (bad[i]) {
      continue;
    }
    const auto rel_err =
        std::abs((calculated_vector_b[i] - current[i]) / (calculated_vector_b[i] + current[i]));
    error = std::max(rel_err, error);
    num_fail = rel_err > tolerance;
  }

  return std::make_pair(error, num_fail);
}

inline double SCALE::compute_final_error(const std::vector<double>& col,
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

inline double SCALE::compute_scale(const std::vector<ThinPixel<double>>& pixels,
                                   const std::vector<double>& weights,
                                   std::size_t offset) noexcept {
  if (pixels.size() == 0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  double sum = 0.0;
  double norm_sum = 0.0;

  for (const auto& p : pixels) {
    const auto w1 = weights[p.bin1_id - offset];
    const auto w2 = weights[p.bin2_id - offset];

    if (!std::isnan(w1) && !std::isnan(w2)) {
      const auto cfx = p.bin1_id != p.bin2_id ? 2.0 : 1.0;
      sum += p.count * cfx;
      norm_sum += (p.count * cfx) / (w1 * w2);
    }
  }

  return std::sqrt(norm_sum / sum);
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
