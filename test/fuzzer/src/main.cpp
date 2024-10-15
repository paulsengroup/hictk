// Copyright (C) 2024 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/common.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <CLI/Error.hpp>
#include <cassert>
#include <cstdio>
#include <exception>
#include <memory>
#include <stdexcept>
#include <tuple>

#include "hictk/fuzzer/cli.hpp"
#include "hictk/fuzzer/config.hpp"
#include "hictk/fuzzer/tools.hpp"
#include "hictk/version.hpp"

using namespace hictk::fuzzer;

extern "C" const char* __asan_default_options() { return "detect_leaks=0"; }  // NOLINT

static void setup_logger_console() {
  auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  //                        [2021-08-12 17:49:34.581] [info]: my log msg
  stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");

  auto main_logger = std::make_shared<spdlog::logger>("main_logger", stderr_sink);

  spdlog::set_default_logger(main_logger);
}

static void setup_logger_console(int verbosity_lvl, bool print_version) {
  spdlog::set_level(spdlog::level::level_enum(verbosity_lvl));
  for (auto& sink : spdlog::default_logger()->sinks()) {
    sink->set_level(spdlog::level::level_enum(verbosity_lvl));
  }

  if (print_version) {
    SPDLOG_INFO(FMT_STRING("[executor] Fuzzing hictk v{}"), hictk::config::version::str());
  }
}

static std::tuple<int, Cli::subcommand, Config> parse_cli_and_setup_logger(Cli& cli) {
  try {
    auto config = cli.parse_arguments();
    const auto subcmd = cli.get_subcommand();
    setup_logger_console(config.verbosity, subcmd == Cli::subcommand::fuzz);
    return std::make_tuple(cli.exit(), subcmd, config);
  } catch (const CLI::ParseError& e) {
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli.exit(e), Cli::subcommand::help, Config());

  } catch (const std::filesystem::filesystem_error& e) {
    SPDLOG_ERROR(FMT_STRING("FAILURE! {}"), e.what());
    return std::make_tuple(1, Cli::subcommand::help, Config());
  } catch (const spdlog::spdlog_ex& e) {
    fmt::print(
        stderr,
        FMT_STRING(
            "FAILURE! An error occurred while setting up the main application logger: {}.\n"),
        e.what());
    return std::make_tuple(1, Cli::subcommand::help, Config());
  }
}

[[nodiscard]] static std::string task_id_to_str(std::uint16_t task_id) {
  if (task_id == 0) {
    return "executor";
  }
  return fmt::to_string(task_id);
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char** argv) noexcept {
  std::unique_ptr<Cli> cli{nullptr};
  auto proc_name = task_id_to_str(0);

  try {
    setup_logger_console();
    cli = std::make_unique<Cli>(argc, argv);
    const auto [ec, subcmd, config] = parse_cli_and_setup_logger(*cli);
    if (ec != 0 || subcmd == Cli::subcommand::help) {
      return ec;
    }

    proc_name = task_id_to_str(config.task_id);
    using sc = Cli::subcommand;
    switch (subcmd) {
      case sc::fuzz:
        return fuzz_subcommand(config);
      case sc::launch_worker:
        return launch_worker_subcommand(config);
      case sc::help:  // NOLINT
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in hictk::fuzzer::main() should be unreachable! "
            "If you see this message, please file an issue on GitHub");
    }
  } catch (const CLI::ParseError& e) {
    assert(cli);
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const std::bad_alloc& err) {
    fmt::print(stderr, FMT_STRING("FAILURE! Unable to allocate enough memory: {}\n"), err.what());
    return 1;
  } catch (const std::exception& e) {
    if (cli) {
      fmt::print(stderr,
                 FMT_STRING("FAILURE! hictk_fuzzer {} [{}] encountered the following error: {}\n"),
                 cli->get_printable_subcommand(), proc_name, e.what());
    } else {
      fmt::print(stderr, FMT_STRING("FAILURE! hictk_fuzzer encountered the following error: {}\n"),
                 e.what());
    }
    return 1;
  } catch (...) {
    fmt::print(
        stderr,
        FMT_STRING("FAILURE! hictk_fuzzer {} [{}] encountered the following error: Caught an "
                   "unhandled exception! "
                   "If you see this message, please file an issue on GitHub.\n"),
        cli->get_printable_subcommand(), proc_name);
    return 1;
  }
  return 0;
}
