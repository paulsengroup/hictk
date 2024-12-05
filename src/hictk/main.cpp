// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/format.h>
#include <fmt/std.h>
#include <spdlog/sinks/callback_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <atomic>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <deque>
#include <memory>
#include <mutex>
#include <string_view>
#include <tuple>
#include <utility>
#include <variant>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"
#include "hictk/tools/tools.hpp"
#include "hictk/version.hpp"

using namespace hictk::tools;

template <std::size_t CAPACITY>
class GlobalLogger {
  //                                              [2021-08-12 17:49:34.581] [info]: my log msg
  static constexpr std::string_view _msg_pattern{"[%Y-%m-%d %T.%e] %^[%l]%$: %v"};
  using HictkLogMsg = std::pair<spdlog::level::level_enum, std::string>;
  std::deque<HictkLogMsg> _msg_buffer{};

  std::mutex _mtx;
  std::atomic<std::size_t> _num_msg_enqueued{};
  std::atomic<bool> _ok{false};

  [[nodiscard]] static std::shared_ptr<spdlog::sinks::stderr_color_sink_mt> init_stderr_sink() {
    auto stderr_sink = std::make_shared<spdlog::sinks::stderr_color_sink_mt>();
    stderr_sink->set_pattern(std::string{_msg_pattern});
    stderr_sink->set_level(spdlog::level::debug);

    return stderr_sink;
  }

  [[nodiscard]] std::shared_ptr<spdlog::sinks::callback_sink_mt> init_callback_sink() {
    if constexpr (CAPACITY != 0 && SPDLOG_ACTIVE_LEVEL <= SPDLOG_LEVEL_WARN) {
      auto callback_sink = std::make_shared<spdlog::sinks::callback_sink_mt>(
          [this](const spdlog::details::log_msg &msg) noexcept { enqueue_msg(msg); });
      callback_sink->set_pattern(std::string{_msg_pattern});
      callback_sink->set_level(spdlog::level::warn);

      return callback_sink;
    }

    return {};
  }

  template <typename... T>
  void print_noexcept(fmt::format_string<T...> fmt, T &&...args) noexcept {
    try {
      fmt::print(stderr, fmt, std::forward<T>(args)...);
    } catch (...) {  // NOLINT
      _ok = false;
    }
  }

  void enqueue_msg(const spdlog::details::log_msg &msg) noexcept {
    if (msg.level < spdlog::level::warn) {
      return;
    }

    ++_num_msg_enqueued;

    try {
      [[maybe_unused]] const std::scoped_lock lck(_mtx);
      if (_msg_buffer.size() == CAPACITY) {
        _msg_buffer.pop_front();
      }
      _msg_buffer.emplace_back(msg.level, std::string{msg.payload.begin(), msg.payload.end()});
    } catch (...) {  // NOLINT
    }
  }

  void replay_warnings() {
    [[maybe_unused]] const std::scoped_lock lck(_mtx);
    if (_msg_buffer.empty()) {
      return;
    }
    auto stderr_sink = init_stderr_sink();
    stderr_sink->set_level(spdlog::level::warn);

    auto logger = std::make_shared<spdlog::logger>("tmp_logger", std::move(stderr_sink));
    logger->set_level(spdlog::level::warn);

    if (_num_msg_enqueued <= _msg_buffer.size()) {
      logger->warn(FMT_STRING("replaying the last {} warning message(s)"),
                   _num_msg_enqueued.load());
    } else {
      logger->warn(FMT_STRING("replaying the last {}/{} warning messages"), _msg_buffer.size(),
                   _num_msg_enqueued.load());
    }
    for (const auto &msg : _msg_buffer) {
      logger->log(msg.first, msg.second);
    }
    _msg_buffer.clear();
  }

  static void reset_logger() noexcept {
    try {
      spdlog::set_default_logger(
          std::make_shared<spdlog::logger>("main_logger", init_stderr_sink()));
    } catch (...) {  // NOLINT
    }
  }

 public:
  GlobalLogger() noexcept {
    try {
      spdlog::set_default_logger(std::make_shared<spdlog::logger>(
          "main_logger", spdlog::sinks_init_list{init_stderr_sink(), init_callback_sink()}));
      _ok = true;
    } catch (const std::exception &e) {
      print_noexcept(FMT_STRING("FAILURE! Failed to setup hictk's logger: {}"), e.what());
    } catch (...) {
      print_noexcept(FMT_STRING("FAILURE! Failed to setup hictk's logger: unknown error"));
    }
  }

  GlobalLogger(const GlobalLogger &other) = delete;
  GlobalLogger(GlobalLogger &&other) noexcept
      : _msg_buffer(std::move(other._msg_buffer)),
        _num_msg_enqueued(other._num_msg_enqueued.load()),
        _ok(other._ok.load()) {
    other._num_msg_enqueued = 0;
    other._ok = false;
  }

  GlobalLogger &operator=(const GlobalLogger &other) = delete;
  GlobalLogger &operator=(GlobalLogger &&other) noexcept {
    if (this == &other) {
      return *this;
    }

    [[maybe_unused]] const auto lck = std::scoped_lock(other._mtx);
    _msg_buffer = std::move(other._msg_buffer);
    _num_msg_enqueued = other._num_msg_enqueued.load();
    _ok = other._ok.load();

    other._num_msg_enqueued = 0;
    other._ok = false;

    return *this;
  }

  ~GlobalLogger() noexcept {
    if (!_ok) {
      reset_logger();
      return;
    }

    try {
      replay_warnings();
    } catch (const std::exception &e) {
      print_noexcept(FMT_STRING("FAILURE! Failed to replay hictk warnings: {}"), e.what());
    } catch (...) {
      print_noexcept(FMT_STRING("FAILURE! Failed to replay hictk warnings: unknown error"));
    }
    reset_logger();
  }

  static void set_level(int lvl) {
    if (auto logger = spdlog::default_logger(); logger) {
      for (auto &sink : logger->sinks()) {
        sink->set_level(std::max(sink->level(), spdlog::level::level_enum{lvl}));
      }
      logger->set_level(spdlog::level::level_enum{lvl});
    }
  }

  void print_welcome_msg() {
    if (_ok) {
      SPDLOG_INFO(FMT_STRING("Running hictk v{}"), hictk::config::version::str());
    }
  }

  [[nodiscard]] constexpr bool ok() const noexcept { return _ok; }

  void clear() noexcept {
    if (_ok) {
      _msg_buffer.clear();
      _num_msg_enqueued = 0;
    }
  }
};

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
              if (subcmd != Cli::subcommand::help && subcmd != Cli::subcommand::dump) {
                logger.print_welcome_msg();
              }
            }
          }
        },
        config);

    return std::make_tuple(cli.exit(), subcmd, config);
  } catch (const CLI::ParseError &e) {
    //  This takes care of formatting and printing error messages (if any)
    return std::make_tuple(cli.exit(e), Cli::subcommand::help, Config());

  } catch (const std::filesystem::filesystem_error &e) {
    SPDLOG_ERROR(FMT_STRING("FAILURE! {}"), e.what());
    return std::make_tuple(1, Cli::subcommand::help, Config());
  } catch (const spdlog::spdlog_ex &e) {
    fmt::print(
        stderr,
        FMT_STRING(
            "FAILURE! An error occurred while setting up the main application logger: {}.\n"),
        e.what());
    return std::make_tuple(1, Cli::subcommand::help, Config());
  }
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
    if (ec != 0 || subcmd == Cli::subcommand::help) {
      local_logger->clear();
      return ec;
    }

    cli->log_warnings();

    using sc = Cli::subcommand;
    switch (subcmd) {
      case sc::balance:
        if (std::holds_alternative<BalanceICEConfig>(config)) {
          return balance_subcmd(std::get<BalanceICEConfig>(config));
        }
        if (std::holds_alternative<BalanceSCALEConfig>(config)) {
          return balance_subcmd(std::get<BalanceSCALEConfig>(config));
        }
        assert(std::holds_alternative<BalanceVCConfig>(config));
        return balance_subcmd(std::get<BalanceVCConfig>(config));
      case sc::convert:
        return convert_subcmd(std::get<ConvertConfig>(config));
      case sc::dump:
        return dump_subcmd(std::get<DumpConfig>(config));
      case sc::fix_mcool:
        return fix_mcool_subcmd(std::get<FixMcoolConfig>(config));
      case sc::load:
        return load_subcmd(std::get<LoadConfig>(config));
      case sc::merge:
        return merge_subcmd(std::get<MergeConfig>(config));
      case sc::metadata:
        return metadata_subcmd(std::get<MetadataConfig>(config));
      case sc::rename_chromosomes:
        return rename_chromosomes_subcmd(std::get<RenameChromosomesConfig>(config));
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
  } catch (const CLI::ParseError &e) {
    assert(cli);
    return cli->exit(e);  //  This takes care of formatting and printing error messages (if any)
  } catch (const std::bad_alloc &err) {
    SPDLOG_CRITICAL(FMT_STRING("FAILURE! Unable to allocate enough memory: {}"), err.what());
    return 1;
  } catch (const spdlog::spdlog_ex &e) {
    fmt::print(stderr,
               FMT_STRING("FAILURE! hictk encountered the following error while logging: {}\n"),
               e.what());
    return 1;
  } catch (const std::exception &e) {
    if (cli) {
      SPDLOG_CRITICAL(FMT_STRING("FAILURE! hictk {} encountered the following error: {}"),
                      cli->get_printable_subcommand(), e.what());
    } else {
      SPDLOG_CRITICAL(FMT_STRING("FAILURE! hictk encountered the following error: {}"), e.what());
    }
    return 1;
  } catch (...) {
    SPDLOG_CRITICAL(FMT_STRING("FAILURE! hictk {} encountered the following error: Caught an "
                               "unhandled exception! "
                               "If you see this message, please file an issue on GitHub."),
                    cli->get_printable_subcommand());
    return 1;
  }
  return 0;
}
