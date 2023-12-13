// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <algorithm>
#include <cstdint>
#include <string>
#include <string_view>
#include <utility>

#include "hictk/file.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

void print(const Pixel<double>& pixel);
void print(const ThinPixel<double>& pixel);

void dump_bins(const File& f, std::string_view range);
void dump_cells(std::string_view uri, std::string_view format);
void dump_chroms(std::string_view uri, std::string_view format, std::uint32_t resolution);
void dump_normalizations(std::string_view uri, std::string_view format, std::uint32_t resolution);
void dump_resolutions(std::string_view uri, std::string_view format, std::uint32_t resolution);

[[nodiscard]] std::pair<std::string, std::string> parse_bedpe(std::string_view line);

template <typename PixelIt>
inline void print_pixels(PixelIt first, PixelIt last) {
  std::for_each(first, last, [&](const auto& pixel) { print(pixel); });
}

}  // namespace hictk::tools
