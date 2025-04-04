// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <stdexcept>
#include <variant>

#include "./config.hpp"

namespace hictk::tools {

[[nodiscard]] inline int run_subcmd([[maybe_unused]] const std::monostate& c) {
  throw std::logic_error(
      "hictk::tools::run_subcmd(const std::monostate& c) should never be called!");
}
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
