// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include "hictk/binary_buffer.hpp"
#include "hictk/pixel.hpp"

namespace hictk::hic::internal {

template <typename Bin1Type, typename Bin2Type, typename CountType>
void HiCBlockReader::read_type1_block(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                      BinaryBuffer &src,
                                      std::vector<ThinPixel<float>> &dest) noexcept {
  using i16 = std::int16_t;
  using i32 = std::int32_t;
  using f32 = float;
  static_assert(std::is_same_v<i16, Bin1Type> || std::is_same_v<i32, Bin1Type>);
  static_assert(std::is_same_v<i16, Bin2Type> || std::is_same_v<i32, Bin2Type>);
  static_assert(std::is_same_v<i16, CountType> || std::is_same_v<f32, CountType>);

  constexpr auto expectedOffsetV7 = (3 * sizeof(i32)) + (2 * sizeof(char));
  constexpr auto expectedOffsetV8plus = expectedOffsetV7 + (2 * sizeof(char));
  std::ignore = expectedOffsetV7;
  std::ignore = expectedOffsetV8plus;
  assert(src() == expectedOffsetV7 || src() == expectedOffsetV8plus);

  const auto expectedNumRecords = dest.size();
  dest.clear();
  const auto numRows = static_cast<i32>(src.read<Bin2Type>());
  for (i32 i = 0; i < numRows; ++i) {
    const auto bin2 = bin2Offset + static_cast<i32>(src.read<Bin2Type>());

    const auto numCols = static_cast<i32>(src.read<Bin1Type>());
    for (i32 j = 0; j < numCols; ++j) {
      const auto bin1 = bin1Offset + static_cast<i32>(src.read<Bin1Type>());

      const auto counts = static_cast<f32>(src.read<CountType>());
      dest.push_back(ThinPixel<float>{static_cast<std::uint64_t>(bin1),
                                      static_cast<std::uint64_t>(bin2), counts});
    }
  }

  std::ignore = expectedNumRecords;
  assert(expectedNumRecords == dest.size());
}

template <typename CountType>
void HiCBlockReader::read_type2_block(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                      BinaryBuffer &src,
                                      std::vector<ThinPixel<float>> &dest) noexcept {
  using i16 = std::int16_t;
  using i32 = std::int32_t;
  using f32 = float;
  static_assert(std::is_same_v<i16, CountType> || std::is_same_v<f32, CountType>);

  const auto nPts = src.read<i32>();
  const auto w = static_cast<i32>(src.read<i16>());

  constexpr auto i16Sentinel = (std::numeric_limits<i16>::lowest)();
  constexpr auto i16Counts = std::is_same_v<i16, CountType>;

  auto isValid = [&](CountType n) {
    return (i16Counts && static_cast<i16>(n) != i16Sentinel) ||
           (!i16Counts && !std::isnan(static_cast<f32>(n)));
  };

  dest.reserve(static_cast<std::size_t>(nPts));
  dest.clear();
  for (i32 i = 0; i < nPts; ++i) {
    const auto count = src.read<CountType>();
    if (!isValid(count)) {
      continue;
    }
    const auto row = i / w;
    const auto col = i - (row * w);
    const auto bin1 = static_cast<std::uint64_t>(bin1Offset) + static_cast<std::uint64_t>(col);
    const auto bin2 = static_cast<std::uint64_t>(bin2Offset) + static_cast<std::uint64_t>(row);

    dest.emplace_back(ThinPixel<float>{bin1, bin2, static_cast<f32>(count)});
  }
}

}  // namespace hictk::hic::internal
