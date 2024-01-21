// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

// IWYU pragma: private, include "hictk/hic.hpp"

#include <parallel_hashmap/btree.h>
#include <parallel_hashmap/phmap.h>

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/chromosome.hpp"
#include "hictk/pixel.hpp"
#include "hictk/reference.hpp"

namespace hictk::hic::internal {

class ExpectedValuesAggregator {
  std::shared_ptr<const BinTable> _bins{};
  std::size_t _num_bins_gw{};

  using CisKey = Chromosome;
  using TransKey = std::pair<CisKey, CisKey>;
  phmap::flat_hash_map<CisKey, double> _cis_sum{};
  phmap::flat_hash_map<TransKey, double> _trans_sum{};

  std::vector<double> _possible_distances{};
  std::vector<double> _actual_distances{};

  std::vector<double> _weights{};
  phmap::btree_map<Chromosome, double> _scaling_factors{};

 public:
  ExpectedValuesAggregator() = default;
  explicit ExpectedValuesAggregator(std::shared_ptr<const BinTable> bins);
  void add(const ThinPixel<float>& p);
  void add(const Pixel<float>& p);

  void compute_density();

  [[nodiscard]] const std::vector<double>& weights() const noexcept;

  [[nodiscard]] double scaling_factor(const Chromosome& chrom) const;
  [[nodiscard]] const phmap::btree_map<Chromosome, double>& scaling_factors() const noexcept;

 private:
  [[nodiscard]] const Reference& chromosomes() const noexcept;

  inline void init_possible_distances();
  void compute_density_cis();
  void compute_density_trans();

  [[nodiscard]] double at(const Chromosome& chrom) const;
  [[nodiscard]] double at(const Chromosome& chrom1, const Chromosome& chrom2) const;

  [[nodiscard]] double& at(const Chromosome& chrom);
  [[nodiscard]] double& at(const Chromosome& chrom1, const Chromosome& chrom2);
};

}  // namespace hictk::hic::internal

#include "./impl/expected_values_aggregator_impl.hpp"  // NOLINT
