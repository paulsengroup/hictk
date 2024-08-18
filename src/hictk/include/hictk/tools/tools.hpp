// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "./config.hpp"

namespace hictk::tools {

[[nodiscard]] int balance_subcmd(const BalanceICEConfig& c);
[[nodiscard]] int balance_subcmd(const BalanceSCALEConfig& c);
[[nodiscard]] int balance_subcmd(const BalanceVCConfig& c);
[[nodiscard]] int convert_subcmd(const ConvertConfig& c);
[[nodiscard]] int dump_subcmd(const DumpConfig& c);
[[nodiscard]] int fix_mcool_subcmd(const FixMcoolConfig& c);
[[nodiscard]] int load_subcmd(const LoadConfig& c);
[[nodiscard]] int merge_subcmd(const MergeConfig& c);
[[nodiscard]] int metadata_subcmd(const MetadataConfig& c);
[[nodiscard]] int rename_chromosomes_subcmd(const RenameChromosomesConfig& c);
[[nodiscard]] int validate_subcmd(const ValidateConfig& c);
[[nodiscard]] int zoomify_subcmd(const ZoomifyConfig& c);

}  // namespace hictk::tools
