// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <spdlog/common.h>
#include <spdlog/logger.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <cassert>
#include <exception>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string_view>

#include "hictk/common.hpp"
#include "hictk/tools/cli.hpp"
#include "hictk/tools/tools.hpp"
#include "hictk/version.hpp"

using namespace hictk::tools;

static void setup_logger_console() {
  auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
  //                        [2021-08-12 17:49:34.581] [info]: my log msg
  stderr_sink->set_pattern("[%Y-%m-%d %T.%e] %^[%l]%$: %v");

  auto main_logger = std::make_shared<spdlog::logger>("main_logger", stderr_sink);

  spdlog::set_default_logger(main_logger);
}

static void setup_logger_console(int verbosity_lvl, bool print_version) {
  for (auto& sink : spdlog::default_logger()->sinks()) {
    sink->set_level(spdlog::level::level_enum(verbosity_lvl));
  }

  if (print_version) {
    spdlog::info(FMT_STRING("Running hictk v{}"), hictk::config::version::str());
  }
}

static std::tuple<int, Cli::subcommand, Config> parse_cli_and_setup_logger(int argc, char** argv) {
  std::unique_ptr<Cli> cli{nullptr};
  try {
    cli = std::make_unique<Cli>(argc, argv);
    auto config = cli->parse_arguments();
    const auto subcmd = cli->get_subcommand();
    std::visit(
        [&](const auto& config_) {
          using T = hictk::remove_cvref_t<decltype(config_)>;
          constexpr auto is_monostate = std::is_same_v<T, std::monostate>;
          if constexpr (!is_monostate) {
            setup_logger_console(config_.verbosity, subcmd != Cli::subcommand::help &&
                                                        subcmd != Cli::subcommand::dump);
          }
        },
        config);

    return std::make_tuple(0, subcmd, config);
  } catch (const CLI::ParseError& e) {
    assert(cli);
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli->exit(e), Cli::subcommand::help, Config());

  } catch (const std::filesystem::filesystem_error& e) {
    spdlog::error(FMT_STRING("FAILURE! {}"), e.what());
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

template <typename... Args>
static void try_log_fatal_error(fmt::format_string<Args...> fmt, Args&&... args) {
  if (spdlog::default_logger()) {
    spdlog::error(fmt, std::forward<Args>(args)...);
  } else {
    fmt::print(stderr, fmt, std::forward<Args>(args)...);
  }
}

// NOLINTNEXTLINE(bugprone-exception-escape)
int main(int argc, char** argv) noexcept {
  std::unique_ptr<Cli> cli{nullptr};
  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  try {
    setup_logger_console();
    const auto [ec, subcmd, config] = parse_cli_and_setup_logger(argc, argv);
    if (ec != 0 || subcmd == Cli::subcommand::help) {
      return ec;
    }
    using sc = Cli::subcommand;
    switch (subcmd) {
      case sc::convert:
        return convert_subcmd(std::get<ConvertConfig>(config));
      case sc::dump:
        return dump_subcmd(std::get<DumpConfig>(config));
      case sc::load:
        return load_subcmd(std::get<LoadConfig>(config));
      case sc::merge:  // NOLINT
        // merge_subcmd(std::get<MergeConfig>(config));
        break;
      case sc::validate:
        return validate_subcmd(std::get<ValidateConfig>(config));
      case sc::zoomify:
        return zoomify_subcmd(std::get<ZoomifyConfig>(config));
      case sc::help:  // NOLINT
        break;
      default:
        throw std::runtime_error(
            "Default branch in switch statement in hictk::main() should be unreachable! "
            "If you see this message, please file an issue on GitHub");
    }
  } catch (const CLI::ParseError& e) {
    assert(cli);
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const std::bad_alloc& err) {
    fmt::print(stderr, FMT_STRING("FAILURE! Unable to allocate enough memory: {}"), err.what());
    return 1;
  } catch (const std::exception& e) {
    if (cli) {
      fmt::print(stderr, FMT_STRING("FAILURE! hictk {} encountered the following error: {}."),
                 cli->get_printable_subcommand(), e.what());
    } else {
      fmt::print(stderr, FMT_STRING("FAILURE! hictk encountered the following error: {}."),
                 e.what());
    }
    return 1;
  } catch (...) {
    fmt::print(stderr,
               FMT_STRING("FAILURE! hictk {} encountered the following error: Caught an "
                          "unhandled exception! "
                          "If you see this message, please file an issue on GitHub."),
               cli->get_printable_subcommand());
    return 1;
  }
  return 0;
}
