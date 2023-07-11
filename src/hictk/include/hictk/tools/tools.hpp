// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "config.hpp"

namespace hictk::tools {

void convert_subcmd(const ConvertConfig& c);
void dump_subcmd(const DumpConfig& c);
void load_subcmd(const LoadConfig& c);
void merge_subcmd(const MergeConfig& c);
[[nodiscard]] int validate_subcmd(const ValidateConfig& c);

}  // namespace hictk::tools
