// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

namespace hictk::hic::internal {

constexpr HiCHeader::operator bool() const noexcept { return footerPosition >= 0; }

}  // namespace hictk::hic::internal
