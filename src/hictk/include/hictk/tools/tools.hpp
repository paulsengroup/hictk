// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "config.hpp"

namespace hictk::tools {

[[nodiscard]] int convert_subcmd(const ConvertConfig& c);
[[nodiscard]] int dump_subcmd(const DumpConfig& c);
[[nodiscard]] int load_subcmd(const LoadConfig& c);
[[nodiscard]] int merge_subcmd(const MergeConfig& c);
[[nodiscard]] int zoomify_subcmd(const ZoomifyConfig& c);

}  // namespace hictk::tools
