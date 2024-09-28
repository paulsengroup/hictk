// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "./common.hpp"

#include <cassert>
#include <cstdint>
#include <filesystem>
#include <utility>
#include <variant>

#include "./init_bin_table.hpp"
#include "./pixel_parser.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/reference.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::tools {

Stats& Stats::operator+=(const Stats& other) {
  std::visit(
      [&](auto& sum_) {
        using T = remove_cvref_t<decltype(sum_)>;

        sum_ += std::get<T>(other.sum);
      },
      sum);
  nnz += other.nnz;

  return *this;
}

[[nodiscard]] PixelParser init_pixel_parser(Format format,
                                            const std::filesystem::path& path_to_interactions,
                                            const std::filesystem::path& path_to_chrom_sizes,
                                            const std::filesystem::path& path_to_bins,
                                            std::uint32_t resolution, std::string_view assembly) {
  assert(format == Format::_4DN || !path_to_chrom_sizes.empty() || !path_to_bins.empty());
  BinTable bins{};

  if (!path_to_bins.empty()) {
    bins = init_bin_table(path_to_chrom_sizes, path_to_bins, resolution);
  } else if (!path_to_chrom_sizes.empty()) {
    bins = BinTable{Reference::from_chrom_sizes(path_to_chrom_sizes), resolution};
  }

  switch (format) {
    case Format::_4DN: {
      if (bins.empty()) {
        assert(resolution != 0);
        return {path_to_interactions, resolution, format, assembly};
      }
      return {path_to_interactions, std::move(bins), format, assembly};
    }
    default:
      return {path_to_interactions, std::move(bins), format, assembly};
  }
  unreachable_code();
}

Format format_from_string(std::string_view s) {
  if (s == "coo") {
    return Format::COO;
  }
  if (s == "bg2") {
    return Format::BG2;
  }
  if (s == "validpairs") {
    return Format::VP;
  }
  assert(s == "4dn");
  return Format::_4DN;
}

}  // namespace hictk::tools
