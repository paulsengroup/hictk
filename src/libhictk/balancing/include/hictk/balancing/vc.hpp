// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <vector>

namespace hictk::balancing {

class VC {
  std::vector<std::uint64_t> _chrom_offsets{};
  std::vector<double> _biases{};
  std::vector<double> _scale{};

  struct Result {
    std::vector<std::uint64_t> offsets{};
    std::vector<double> scales{};
    std::vector<double> weights{};
  };

 public:
  enum class Type : std::uint_fast8_t { cis, trans, gw };

  struct Params {};

  template <typename File>
  explicit VC(const File& f, Type type = Type::gw, const Params& params = {});
  template <typename PixelIt>
  VC(PixelIt first, PixelIt last, const BinTable& bins, const Params& params = {});

  [[nodiscard]] Weights get_weights(bool rescale = true) const;
  [[nodiscard]] const std::vector<double>& get_scale() const noexcept;

 private:
  template <typename File>
  [[nodiscard]] static auto compute_cis(const File& f) -> Result;
  template <typename File>
  [[nodiscard]] static auto compute_trans(const File& f) -> Result;
  template <typename File>
  [[nodiscard]] static auto compute_gw(const File& f) -> Result;
};

}  // namespace hictk::balancing

#include "./impl/vc_impl.hpp"  //NOLINT
