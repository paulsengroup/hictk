// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace hictk::balancing::internal {

constexpr bool VectorOfAtomicDecimals::overflows(double n) const noexcept { return n > _max_value; }

}  // namespace hictk::balancing::internal
