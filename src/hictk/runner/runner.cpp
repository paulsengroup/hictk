// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/runner.hpp"

#include <stdexcept>
#include <variant>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/telemetry.hpp"
#include "include/hictk/tools/tools.hpp"

namespace hictk::tools {

int run_subcommand(Cli::subcommand subcmd, const Config &config) {
  if (subcmd == Cli::subcommand::none) {
    throw std::runtime_error(
        "run_subcommand() was called with subcommand::none: this should never happen! "
        "If you see this message, please file an issue on GitHub");
  }

  auto *tracer = Tracer::instance();
  return std::visit(
      [&](const auto &c) {
        using ScopedSpan = decltype(tracer->get_scoped_span(subcmd, c));
        auto span = !!tracer ? tracer->get_scoped_span(subcmd, c) : ScopedSpan{};
        const auto ec = run_subcmd(c);
        if (ec == 0 && span.has_value()) {
          span->set_status(Tracer::StatusCode::kOk);
        }

        return ec;
      },
      config);
}

void try_tear_down_telemetry_reporter() noexcept { Tracer::tear_down_instance(); }

}  // namespace hictk::tools
