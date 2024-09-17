// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <readerwriterqueue/readerwriterqueue.h>

#include <cassert>
#include <cstdint>
#include <string_view>
#include <variant>
#include <vector>

#include "hictk/bin_table.hpp"
#include "hictk/pixel.hpp"
#include "hictk/type_traits.hpp"

namespace hictk::tools {

template <typename N>
using PixelQueue = moodycamel::BlockingReaderWriterQueue<ThinPixel<N>>;

using PixelQueueVar = std::variant<PixelQueue<std::int32_t>, PixelQueue<float>, PixelQueue<double>>;

using IntBuff = std::vector<ThinPixel<std::int32_t>>;
using FPBuff = std::vector<ThinPixel<double>>;
using PixelBuffer = std::variant<IntBuff, FPBuff>;

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

[[nodiscard]] inline bool line_is_header(std::string_view line) {
  return !line.empty() && line.front() == '#';
}

struct Stats {
  std::variant<std::uint64_t, double> sum{0.0};
  std::uint64_t nnz{};

  inline Stats& operator+=(const Stats& other) {
    std::visit(
        [&](auto& sum_) {
          using T = remove_cvref_t<decltype(sum_)>;

          sum_ += std::get<T>(other.sum);
        },
        sum);
    nnz += other.nnz;

    return *this;
  }

  template <typename N>
  inline Stats(N sum_, std::uint64_t nnz_) : nnz(nnz_) {
    if constexpr (std::is_floating_point_v<N>) {
      sum = static_cast<double>(sum_);
    } else {
      sum = static_cast<std::uint64_t>(sum_);
    }
  }
};

}  // namespace hictk::tools
