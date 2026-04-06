// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "hictk/hic/common.hpp"

namespace hictk::hic {

constexpr auto File::matrix_type() const noexcept -> MatrixType { return _type; }

constexpr auto File::matrix_unit() const noexcept -> MatrixUnit { return _unit; }

}  // namespace hictk::hic
