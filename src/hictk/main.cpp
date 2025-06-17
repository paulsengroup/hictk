// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <CLI/Error.hpp>
#include <exception>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/logger.hpp"
#include "hictk/tools/runner.hpp"
#include "hictk/type_traits.hpp"

using namespace hictk::tools;

// NOLINTNEXTLINE(*-err58-cpp, *-avoid-non-const-global-variables, *-avoid-magic-numbers)
static auto global_logger = std::make_unique<GlobalLogger<256>>();

static auto acquire_global_logger() noexcept { return std::move(global_logger); }

template <typename Logger>
static std::tuple<int, Cli::subcommand, Config> parse_cli_and_setup_logger(Cli &cli,
                                                                           Logger &logger) {
  try {
    auto config = cli.parse_arguments();
    const auto subcmd = cli.get_subcommand();
    std::visit(
        [&](const auto &config_) {
          using T = hictk::remove_cvref_t<decltype(config_)>;
          if constexpr (!std::is_same_v<T, std::monostate>) {
            if (logger.ok()) {
              logger.set_level(config_.verbosity);  // NOLINT
              if (subcmd != Cli::subcommand::none && subcmd != Cli::subcommand::dump) {
                logger.print_welcome_msg();
              }
            }
          }
        },
        config);

    return std::make_tuple(cli.exit(), subcmd, config);
  } catch (const CLI::ParseError &e) {
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli.exit(e), Cli::subcommand::none, Config());

  } catch (const std::filesystem::filesystem_error &e) {
    SPDLOG_ERROR(FMT_STRING("FAILURE! {}"), e.what());
    return std::make_tuple(1, Cli::subcommand::none, Config());
  } catch (const spdlog::spdlog_ex &e) {
    fmt::print(
        stderr,
        FMT_STRING(
            "FAILURE! An error occurred while setting up the main application logger: {}.\n"),
        e.what());
    return std::make_tuple(1, Cli::subcommand::none, Config());
  }
}

[[nodiscard]] static std::string generate_command_name(const std::unique_ptr<Cli> &cli) {
  if (cli) {
    return fmt::format(FMT_STRING("hictk {}"), cli->get_printable_subcommand());
  }

  return "hictk";
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char **argv) noexcept {
  std::unique_ptr<Cli> cli{nullptr};

  try {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    auto local_logger = acquire_global_logger();
    cli = std::make_unique<Cli>(argc, argv);
    const auto [ec, subcmd, config] = parse_cli_and_setup_logger(*cli, *local_logger);
    if (ec != 0 || subcmd == Cli::subcommand::none) {
      local_logger->clear();
      return ec;
    }

    cli->log_warnings();

    const auto ec_ = run_subcommand(subcmd, config);
    try_tear_down_telemetry_reporter();
    return ec_;
  } catch (const CLI::ParseError &e) {
    if (cli) {
      //  This takes care of formatting and printing error messages (if any)
      return cli->exit(e);
    }
    SPDLOG_CRITICAL("FAILURE! An unknown error occurred while parsing CLI arguments.");
  } catch (const std::bad_alloc &e) {
    SPDLOG_CRITICAL(FMT_STRING("FAILURE! Unable to allocate enough memory: {}"), e.what());
  } catch (const spdlog::spdlog_ex &e) {
    fmt::print(stderr,
               FMT_STRING("FAILURE! {} encountered the following error while logging: {}\n"),
               generate_command_name(cli), e.what());
  } catch (const std::exception &e) {
    SPDLOG_CRITICAL(FMT_STRING("FAILURE! {} encountered the following error: {}"),
                    generate_command_name(cli), e.what());
  } catch (...) {
    SPDLOG_CRITICAL(
        FMT_STRING("FAILURE! {} encountered the following error: Caught an unhandled exception! If "
                   "you see this message, please file an issue on GitHub."),
        generate_command_name(cli));
  }
  try_tear_down_telemetry_reporter();
  return 1;
}
