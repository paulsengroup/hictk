// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/compile.h>
#include <fmt/format.h>

#include "hictk/cooler.hpp"

namespace hictk::tools {

template <bool join>
inline void print(const Pixel<std::int64_t>& pixel) {
  if constexpr (join) {
    fmt::print(FMT_COMPILE("{:bg2}\n"), pixel);
  } else {
    fmt::print(FMT_COMPILE("{:raw}\n"), pixel);
  }
}

template <bool join>
inline void print(const Pixel<double>& pixel) {
  if constexpr (join) {
    fmt::print(FMT_COMPILE("{:bg2}\t{:g}\n"), pixel.coords, pixel.count);
  } else {
    fmt::print(FMT_COMPILE("{:raw}\t{:g}\n"), pixel.coords, pixel.count);
  }
}

template <typename File>
inline void dump_chroms(const File& f, std::string_view range) {
  if (range == "all") {
    for (const Chromosome& chrom : f.chromosomes()) {
      fmt::print(FMT_COMPILE("{:s}\t{:d}\n"), chrom.name(), chrom.size());
    }
    return;
  }

  const auto coords = GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range});
  auto it = f.chromosomes().find(coords.chrom());
  if (it != f.chromosomes().end()) {
    fmt::print(FMT_COMPILE("{:s}\t{:d}\n"), it->name(), it->size());
  }
}

template <typename File>
inline void dump_bins(const File& f, std::string_view range) {
  if (range == "all") {
    for (const auto& bin : f.bins()) {
      fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
    }
    return;
  }

  const auto coords = GenomicInterval::parse_ucsc(f.chromosomes(), std::string{range});
  auto [first_bin, last_bin] = f.bins().find_overlap(coords);
  std::for_each(first_bin, last_bin, [](const Bin& bin) {
    fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\n"), bin.chrom().name(), bin.start(), bin.end());
  });
}

[[nodiscard]] inline std::pair<std::string, std::string> parse_bedpe(std::string_view line) {
  auto next_token = [&]() {
    assert(!line.empty());
    const auto pos1 = line.find('\t');
    const auto pos2 = line.find('\t', pos1 + 1);
    const auto pos3 = line.find('\t', pos2 + 1);

    auto tok = std::string{line.substr(0, pos3)};
    tok[pos1] = ':';
    tok[pos2] = '-';
    line.remove_prefix(pos3 + 1);
    return tok;
  };

  return std::make_pair(next_token(), next_token());
}
}  // namespace hictk::tools
