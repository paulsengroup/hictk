// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "./config.hpp"

namespace hictk::tools {

[[nodiscard]] int run_subcmd(const BalanceICEConfig& c);
[[nodiscard]] int run_subcmd(const BalanceSCALEConfig& c);
[[nodiscard]] int run_subcmd(const BalanceVCConfig& c);
[[nodiscard]] int run_subcmd(const ConvertConfig& c);
[[nodiscard]] int run_subcmd(const DumpConfig& c);
[[nodiscard]] int run_subcmd(const FixMcoolConfig& c);
[[nodiscard]] int run_subcmd(const LoadConfig& c);
[[nodiscard]] int run_subcmd(const MergeConfig& c);
[[nodiscard]] int run_subcmd(const MetadataConfig& c);
[[nodiscard]] int run_subcmd(const RenameChromosomesConfig& c);
[[nodiscard]] int run_subcmd(const ValidateConfig& c);
[[nodiscard]] int run_subcmd(const ZoomifyConfig& c);

}  // namespace hictk::tools
