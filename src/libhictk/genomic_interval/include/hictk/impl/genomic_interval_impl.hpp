// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>

namespace hictk {

constexpr std::uint32_t GenomicInterval::start() const noexcept { return _start; }
constexpr std::uint32_t GenomicInterval::end() const noexcept { return _end; }
constexpr std::uint32_t GenomicInterval::size() const noexcept { return _end - _start; }
constexpr bool GenomicInterval::empty() const noexcept { return size() == 0; }

}  // namespace hictk
