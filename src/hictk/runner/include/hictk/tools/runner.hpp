// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {
[[nodiscard]] int run_subcommand(Cli::subcommand subcmd, const Config &config);

void try_tear_down_telemetry_reporter() noexcept;

}  // namespace hictk::tools
