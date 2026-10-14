// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "hictk/bin_table.hpp"
#include "hictk/cooler/cooler.hpp"
#include "hictk/file.hpp"
#include "hictk/hic.hpp"

namespace hictk::balancing::internal {

void check_storage_mode(const cooler::File& clr);
constexpr void check_storage_mode(const hic::File&) noexcept {}
void check_storage_mode(const File& f);

void check_bin_type(const BinTable& bins);

}  // namespace hictk::balancing::internal
