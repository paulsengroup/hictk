// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <string_view>

#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

enum class Format { COO, BG2, VP, _4DN };
[[nodiscard]] inline Format format_from_string(std::string_view s) {
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

template <typename N>
[[nodiscard]] inline ThinPixel<N> parse_pixel(const BinTable& bins, std::string_view line,
                                              Format format) {
  ThinPixel<N> pixel{};
  switch (format) {
    case Format::COO:
      pixel = ThinPixel<N>::from_coo(bins, line);
      break;
    case Format::BG2:
      pixel = Pixel<N>::from_bg2(bins, line).to_thin();
      break;
    case Format::VP:
      pixel = Pixel<N>::from_validpair(bins, line).to_thin();
      break;
    case Format::_4DN:
      pixel = Pixel<N>::from_4dn_pairs(bins, line).to_thin();
      break;
  }
  if (pixel.bin1_id > pixel.bin2_id) {
    std::swap(pixel.bin1_id, pixel.bin2_id);
  }
}

[[nodiscard]] inline bool line_is_header(std::string_view line) {
  return !line.empty() && line.front() == '#';
}

template <typename N>
struct CoolerChunk {
  std::string uri;
  cooler::PixelSelector::iterator<N> first;
  cooler::PixelSelector::iterator<N> last;
};

}  // namespace hictk::tools
