// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>

namespace hictk {

constexpr std::uint64_t Bin::id() const noexcept { return _id; }
constexpr std::uint32_t Bin::rel_id() const noexcept { return _rel_id; }
constexpr std::uint32_t Bin::start() const noexcept { return _interval.start(); }
constexpr std::uint32_t Bin::end() const noexcept { return _interval.end(); }

constexpr bool Bin::has_null_id() const noexcept { return id() == null_id; }

}  // namespace hictk
