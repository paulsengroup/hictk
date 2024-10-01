// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/reference.hpp"

namespace hictk {

class MultiResFile {
  std::string _path{};
  hic::MatrixType _type{};
  hic::MatrixUnit _unit{};
  Reference _chroms{};
  std::vector<std::uint32_t> _resolutions{};
  std::string _format{};
  std::uint8_t _format_version{};
  BinTable::Type _bin_type{};

 public:
  using QUERY_TYPE = hictk::GenomicInterval::Type;

  explicit MultiResFile(const cooler::MultiResFile& mclr);
  explicit MultiResFile(const hic::File& hf);
  explicit MultiResFile(std::string uri, hic::MatrixType type_ = hic::MatrixType::observed,
                        hic::MatrixUnit unit_ = hic::MatrixUnit::BP);

  [[nodiscard]] std::string path() const;

  [[nodiscard]] bool is_hic() const noexcept;
  [[nodiscard]] bool is_mcool() const noexcept;

  [[nodiscard]] constexpr hic::MatrixType matrix_type() const noexcept;
  [[nodiscard]] constexpr hic::MatrixUnit matrix_unit() const noexcept;
  [[nodiscard]] std::string_view format() const noexcept;
  [[nodiscard]] constexpr std::uint8_t version() const noexcept;
  [[nodiscard]] constexpr BinTable::Type bin_type() const noexcept;
  [[nodiscard]] constexpr const std::vector<std::uint32_t>& resolutions() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;

  [[nodiscard]] File open(std::uint32_t resolution) const;
};

}  // namespace hictk

#include "./impl/multires_file_impl.hpp"  // NOLINT
