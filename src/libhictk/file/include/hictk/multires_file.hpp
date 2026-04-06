// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include "hictk/balancing/methods.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/multires_cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"
#include "hictk/reference.hpp"

namespace hictk {

class MultiResFile {
 public:
  using QUERY_TYPE = GenomicInterval::Type;
  enum class Format : std::uint_fast8_t { UNKNOWN, MCOOL, HIC };

 private:
  std::string _path{};
  hic::MatrixType _type{hic::MatrixType::observed};
  hic::MatrixUnit _unit{hic::MatrixUnit::BP};
  Format _format{Format::UNKNOWN};
  Reference _chroms{};
  std::vector<std::uint32_t> _resolutions{};
  mutable std::optional<std::pair<std::string, std::vector<balancing::Method>>> _normalizations{};
  std::uint8_t _format_version{};
  BinTable::Type _bin_type{BinTable::Type::fixed};

 public:
  explicit MultiResFile(const cooler::MultiResFile& mclr);
  explicit MultiResFile(const hic::File& hf);
  explicit MultiResFile(std::string uri, hic::MatrixType type_ = hic::MatrixType::observed,
                        hic::MatrixUnit unit_ = hic::MatrixUnit::BP);

  [[nodiscard]] std::string path() const;

  [[nodiscard]] constexpr bool is_hic() const noexcept;
  [[nodiscard]] constexpr bool is_mcool() const noexcept;

  [[nodiscard]] constexpr hic::MatrixType matrix_type() const noexcept;
  [[nodiscard]] constexpr hic::MatrixUnit matrix_unit() const noexcept;
  [[nodiscard]] constexpr std::string_view format() const noexcept;
  [[nodiscard]] constexpr std::uint8_t version() const noexcept;
  [[nodiscard]] constexpr BinTable::Type bin_type() const noexcept;
  [[nodiscard]] constexpr const std::vector<std::uint32_t>& resolutions() const noexcept;
  [[nodiscard]] const Reference& chromosomes() const noexcept;
  [[nodiscard]] const std::vector<balancing::Method>& avail_normalizations(
      std::string_view policy = "union") const;

  [[nodiscard]] File open(std::uint32_t resolution) const;
};

}  // namespace hictk

#include "./impl/multires_file_impl.hpp"  // NOLINT
