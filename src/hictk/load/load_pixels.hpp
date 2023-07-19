// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include "./common.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/common.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

template <typename N>
inline void read_batch(const BinTable& bins, std::vector<ThinPixel<N>>& buffer, Format format) {
  buffer.clear();
  std::string line{};
  try {
    while (std::getline(std::cin, line)) {
      if (line_is_header(line)) {
        continue;
      }
      buffer.emplace_back(parse_pixel<N>(bins, line, format));
      if (buffer.size() == buffer.capacity()) {
        return;
      }
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("encountered error while processing the following line:\n"
                               "\"{}\"\n"
                               "Cause: {}"),
                    line, e.what()));
  }
}

template <typename N>
inline void ingest_pixels_sorted(cooler::File&& clr, Format format, std::size_t batch_size,
                                 bool validate_pixels) {
  std::vector<ThinPixel<N>> buffer(batch_size);

  std::size_t i = 0;
  try {
    for (; !std::cin.eof(); ++i) {
      spdlog::debug(FMT_STRING("processing chunk #{}..."), i + 1);
      read_batch(clr.bins(), buffer, format);
      clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
      buffer.clear();
    }
    if (!buffer.empty()) {
      clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
    }
  } catch (const std::exception& e) {
    const auto i0 = i * buffer.capacity();
    const auto i1 = i0 + buffer.size();
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while processing chunk {}-{}: {}"), i0, i1, e.what()));
  }
}

template <typename N>
inline std::string ingest_pixels_unsorted(cooler::File&& clr, std::vector<ThinPixel<N>>& buffer,
                                          Format format, bool validate_pixels) {
  assert(buffer.capacity() != 0);

  read_batch(clr.bins(), buffer, format);

  if (buffer.empty()) {
    assert(std::cin.eof());
    return "";
  }

  std::sort(buffer.begin(), buffer.end());
  clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
  buffer.clear();

  return clr.uri();
}
}  // namespace hictk::tools
