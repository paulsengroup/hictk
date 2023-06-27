// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "hictk/tools/config.hpp"

namespace hictk::tools {
void hic_to_cool(const ConvertConfig& c);
void cool_to_hic(const ConvertConfig& c);
}  // namespace hictk::tools
