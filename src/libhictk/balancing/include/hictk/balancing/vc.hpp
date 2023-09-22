// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <variant>
#include <vector>

namespace hictk::balancing {

template <typename N>
class VC {
  std::variant<std::vector<double>, std::vector<std::int64_t>> _rowsum{};
  std::variant<std::int64_t, double> _sum{};
  double _norm_sum{};

 public:
  template <typename PixelIt>
  VC(PixelIt first_pixel, PixelIt last_pixel, std::size_t num_rows, std::size_t binid_offset = 0);

  [[nodiscard]] std::vector<double> get_weights() const;
};
}  // namespace hictk::balancing

#include "./impl/vc_impl.hpp"
