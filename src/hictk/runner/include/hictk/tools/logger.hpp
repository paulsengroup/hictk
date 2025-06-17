// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <spdlog/sinks/callback_sink.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/spdlog.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <deque>
#include <exception>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <utility>

#include "hictk/version.hpp"

namespace hictk::tools {

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

}  // namespace hictk::tools
