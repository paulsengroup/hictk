// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <variant>

namespace coolerpp::internal {

// Variants are listed in order from the most common to the least common for perf. reasons
// clang-format off
using NumericVariant =
    std::variant<
        std::uint32_t,
        std::int32_t,
        double,
        std::uint8_t,
        std::uint16_t,
        std::uint64_t,
        std::int8_t,
        std::int16_t,
        std::int64_t,
        float,
        long double>;
// clang-format on

}  // namespace coolerpp::internal
