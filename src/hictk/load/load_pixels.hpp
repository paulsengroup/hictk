// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "./common.hpp"
#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/hic/file_writer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::tools {

template <typename N>
inline Stats read_batch(const BinTable& bins, std::vector<ThinPixel<N>>& buffer, Format format,
                        std::int64_t offset) {
  buffer.clear();
  Stats stats{N{}, 0};
  std::string line{};
  try {
    while (std::getline(std::cin, line)) {
      if (line_is_header(line)) {
        continue;
      }
      const auto& p = buffer.emplace_back(parse_pixel<N>(bins, line, format, offset));
      stats.nnz++;
      if constexpr (std::is_floating_point_v<N>) {
        std::get<double>(stats.sum) += conditional_static_cast<double>(p.count);
      } else {
        std::get<std::uint64_t>(stats.sum) += conditional_static_cast<std::uint64_t>(p.count);
      }
      if (buffer.size() == buffer.capacity()) {
        return stats;
      }
    }
  } catch (const std::exception& e) {
    throw std::runtime_error(
        fmt::format(FMT_STRING("encountered error while processing the following line:\n"
                               "\"{}\"\n"
                               "Cause: {}"),
                    line, e.what()));
  }

  return stats;
}

template <typename N>
[[nodiscard]] inline Stats ingest_pixels_sorted(
    cooler::File&& clr,  // NOLINT(*-rvalue-reference-param-not-moved)
    Format format, std::int64_t offset, std::size_t batch_size, bool validate_pixels) {
  std::vector<ThinPixel<N>> buffer(batch_size);

  std::size_t i = 0;
  Stats stats{N{}, 0};
  try {
    for (; !std::cin.eof(); ++i) {
      SPDLOG_INFO(FMT_STRING("processing chunk #{}..."), i + 1);
      stats += read_batch(clr.bins(), buffer, format, offset);
      clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
      buffer.clear();
    }
    assert(!buffer.empty());
  } catch (const std::exception& e) {
    const auto i0 = i * buffer.capacity();
    const auto i1 = i0 + buffer.size();
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while processing chunk {}-{}: {}"), i0, i1, e.what()));
  }

  return stats;
}

template <typename N>
[[nodiscard]] inline Stats ingest_pixels_unsorted(
    cooler::File&& clr,  // NOLINT(*-rvalue-reference-param-not-moved)
    std::vector<ThinPixel<N>>& buffer, Format format, std::int64_t offset, bool validate_pixels) {
  assert(buffer.capacity() != 0);

  auto stats = read_batch(clr.bins(), buffer, format, offset);

  if (buffer.empty()) {
    assert(std::cin.eof());
    return {N{}, 0};
  }

  std::sort(buffer.begin(), buffer.end());
  clr.append_pixels(buffer.begin(), buffer.end(), validate_pixels);
  buffer.clear();

  clr.flush();
  return stats;
}

[[nodiscard]] inline Stats ingest_pixels(
    hic::internal::HiCFileWriter&& hf,  // NOLINT(*-rvalue-reference-param-not-moved)
    std::vector<ThinPixel<float>>& buffer, Format format, std::int64_t offset) {
  assert(buffer.capacity() != 0);

  std::size_t i = 0;
  Stats stats{0.0, 0};
  try {
    auto t0 = std::chrono::steady_clock::now();
    const auto& bins = hf.bins(hf.resolutions().front());
    for (; !std::cin.eof(); ++i) {
      stats += read_batch(bins, buffer, format, offset);

      if (buffer.empty()) {
        assert(std::cin.eof());
        break;
      }

      const auto t1 = std::chrono::steady_clock::now();
      const auto delta =
          static_cast<double>(
              std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count()) /
          1000.0;
      t0 = t1;
      SPDLOG_INFO(FMT_STRING("preprocessing chunk #{} at {:.0f} pixels/s..."), i + 1,
                  double(buffer.size()) / delta);
      hf.add_pixels(bins.bin_size(), buffer.begin(), buffer.end());
      buffer.clear();
    }
    hf.serialize();
    assert(buffer.empty());
    return stats;
  } catch (const std::exception& e) {
    const auto i0 = i * buffer.capacity();
    const auto i1 = i0 + buffer.size();
    throw std::runtime_error(fmt::format(
        FMT_STRING("an error occurred while processing chunk {}-{}: {}"), i0, i1, e.what()));
  }
}

}  // namespace hictk::tools
