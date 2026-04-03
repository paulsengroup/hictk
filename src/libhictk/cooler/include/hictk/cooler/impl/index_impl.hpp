// Copyright (C) 2022 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>

namespace hictk::cooler {

constexpr std::uint64_t Index::nnz() const noexcept { return _nnz; }

}  // namespace hictk::cooler
