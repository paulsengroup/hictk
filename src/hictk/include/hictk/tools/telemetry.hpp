// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <spdlog/spdlog.h>

#include <cstdint>
#include <optional>

#include "hictk/tools/cli.hpp"

#ifdef HICTK_ENABLE_TELEMETRY

#ifdef OTEL_INTERNAL_LOG_LEVEL
#undef OTEL_INTERNAL_LOG_LEVEL
#endif

#if SPDLOG_ACTIVE_LEVEL < SPDLOG_LEVEL_INFO
#define OTEL_INTERNAL_LOG_LEVEL OTEL_INTERNAL_LOG_DEBUG
#else
#define OTEL_INTERNAL_LOG_LEVEL OTEL_INTERNAL_LOG_NONE
#endif

#include <fmt/format.h>
#include <opentelemetry/exporters/ostream/span_exporter_factory.h>
#include <opentelemetry/exporters/otlp/otlp_http_exporter_factory.h>
#include <opentelemetry/exporters/otlp/otlp_http_exporter_options.h>
#include <opentelemetry/sdk/common/global_log_handler.h>
#include <opentelemetry/sdk/metrics/aggregation/default_aggregation.h>
#include <opentelemetry/sdk/metrics/export/periodic_exporting_metric_reader.h>
#include <opentelemetry/sdk/metrics/meter_provider.h>
#include <opentelemetry/sdk/trace/batch_span_processor_factory.h>
#include <opentelemetry/sdk/trace/batch_span_processor_options.h>
#include <opentelemetry/sdk/trace/simple_processor_factory.h>
#include <opentelemetry/sdk/trace/tracer_provider.h>
#include <opentelemetry/sdk/trace/tracer_provider_factory.h>
#include <opentelemetry/trace/provider.h>
#include <parallel_hashmap/btree.h>

#include <cassert>
#include <chrono>
#include <ciso646>
#include <cstdlib>
#include <exception>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/tools/config.hpp"

namespace hictk::tools {

template <typename SpanPtr, typename ScopePtr>
struct ScopedSpan {
  SpanPtr span;
  ScopePtr scope;

  ScopedSpan(SpanPtr span_, ScopePtr scope_) noexcept
      : span(std::move(span_)), scope(std::move(scope_)) {}

  template <typename Config>
  void update_attributes(const Config& c) noexcept {
    static_assert(!std::is_same_v<Config, std::monostate>);
    try {
      if (span) {
        update_tracer_attributes(span, c);
      }
    } catch (...) {  // NOLINT
    }
  }
  void set_status(opentelemetry::trace::StatusCode s) noexcept {
    try {
      if (span) {
        span->SetStatus(s);
      }

    } catch (...) {  // NOLINT
    }
  }

 private:
  static void update_tracer_attributes(SpanPtr& span, const BalanceICEConfig& c) {
    span->SetAttribute("input-format", infer_input_format(c.path_to_input));
  }

  static void update_tracer_attributes(SpanPtr& span, const BalanceSCALEConfig& c) {
    span->SetAttribute("input-format", infer_input_format(c.path_to_input));
  }

  static void update_tracer_attributes(SpanPtr& span, const BalanceVCConfig& c) {
    span->SetAttribute("input-format", infer_input_format(c.path_to_input));
  }

  static void update_tracer_attributes(SpanPtr& span, const ConvertConfig& c) {
    span->SetAttribute("input-format", c.input_format);
    span->SetAttribute("output-format", c.output_format);
  }

  static void update_tracer_attributes(SpanPtr& span, const DumpConfig& c) {
    auto input_format = infer_input_format(c.uri);
    if (input_format == "mcool" && c.resolution.has_value()) {
      input_format = "cool";
    }
    span->SetAttribute("input-format", input_format);
    span->SetAttribute("table", c.table);
  }

  static void update_tracer_attributes(SpanPtr& span, const LoadConfig& c) {
    span->SetAttribute("output-format", c.output_format);
  }

  static void update_tracer_attributes(SpanPtr& span, const MergeConfig& c) {
    phmap::btree_set<std::string> input_formats{};
    for (const auto& f : c.input_files) {
      input_formats.emplace(infer_input_format(f));
    }
    span->SetAttribute("input-formats",
                       fmt::format(FMT_STRING("{}"), fmt::join(input_formats, ",")));
    span->SetAttribute("output-format", c.output_format);
  }

  static void update_tracer_attributes(SpanPtr& span, const MetadataConfig& c) {
    span->SetAttribute("input-format", c.input_format);
    span->SetAttribute("output-format", c.output_format);
  }

  static void update_tracer_attributes([[maybe_unused]] SpanPtr& span,
                                       [[maybe_unused]] const FixMcoolConfig& c) noexcept {}

  static void update_tracer_attributes(SpanPtr& span, const RenameChromosomesConfig& c) {
    span->SetAttribute("input-format", infer_input_format(c.uri));
  }

  static void update_tracer_attributes(SpanPtr& span, const ValidateConfig& c) {
    try {
      span->SetAttribute("input-format", infer_input_format(c.uri));
    } catch (...) {
      span->SetAttribute("input-format", "unknown");
    }
    span->SetAttribute("output-format", c.output_format);
  }

  static void update_tracer_attributes(SpanPtr& span, const ZoomifyConfig& c) {
    auto input_format = infer_input_format(c.path_to_input);
    if (input_format == "mcool" && c.resolutions.size() == 1) {
      input_format = "cool";
    }
    span->SetAttribute("input-format", input_format);
  }
};

class OpenTelemetryLogHandler : public opentelemetry::sdk::common::internal_log::LogHandler {
  void Handle([[maybe_unused]] opentelemetry::sdk::common::internal_log::LogLevel level,
              [[maybe_unused]] const char* file, [[maybe_unused]] int line, const char* msg,
              [[maybe_unused]] const opentelemetry::sdk::common::AttributeMap& attributes) noexcept
      override {
    if (msg) {
      SPDLOG_DEBUG(FMT_STRING("opentelemetry [{}]: {}"),
                   opentelemetry::sdk::common::internal_log::LevelToString(level), msg);
    }
  }
};

class Tracer {
  using ProviderPtr = decltype(opentelemetry::trace::Provider::GetTracerProvider());
  using TracerPtr = decltype(std::declval<ProviderPtr>()->GetTracer(""));

  using SpanT = remove_cvref_t<decltype(std::declval<TracerPtr>()->StartSpan(""))>;
  using ScopeT = remove_cvref_t<decltype(opentelemetry::trace::Tracer::WithActiveSpan(
      std::declval<std::add_lvalue_reference_t<SpanT>>()))>;

  ProviderPtr _provider{};
  TracerPtr _tracer{};

 public:
  using StatusCode = opentelemetry::trace::StatusCode;

  // NOLINTNEXTLINE(*-mt-unsafe)
  Tracer() noexcept {
    // NOLINTNEXTLINE(*-implicit-bool-conversion, *-mt-unsafe)
    if (!!std::getenv("HICTK_NO_TELEMETRY")) {
      SPDLOG_DEBUG(
          "HICTK_NO_TELEMETRY found in environment variable list: no telemetry information will be "
          "collected.");
      return;
    }

    if (init_remote_telemetry_tracer()) {
      try {
        init_opentelemetry_logger_once();
        _provider = opentelemetry::trace::Provider::GetTracerProvider();
        _tracer = _provider->GetTracer("hictk");
        register_reset_opentelemetry_tracer_at_exit();
      } catch (...) {  // NOLINT
      }
    }
  }

  Tracer(const Tracer&) = delete;
  Tracer(Tracer&&) noexcept = default;

  ~Tracer() noexcept {
    if (!!_tracer) {
      // NOLINTNEXTLINE(*-avoid-magic-numbers)
      _tracer->ForceFlush(std::chrono::milliseconds{500});
    }
  }

  Tracer& operator=(const Tracer&) = delete;
  Tracer& operator=(Tracer&&) noexcept = default;

  template <typename Config>
  [[nodiscard]] auto get_scoped_span(Cli::subcommand subcmd, const Config& config,
                                     StatusCode default_status_code = StatusCode::kError) noexcept
      -> std::optional<ScopedSpan<SpanT, ScopeT>> {
    if (!_tracer) {
      return {};
    }

    std::string subcmd_str{Cli::subcommand_to_str(subcmd)};
    if constexpr (std::is_same_v<BalanceICEConfig, Config>) {
      subcmd_str += "-ice";
    } else if constexpr (std::is_same_v<BalanceSCALEConfig, Config>) {
      subcmd_str += "-scale";
    } else if constexpr (std::is_same_v<BalanceVCConfig, Config>) {
      subcmd_str += "-vc";
    }

    auto span = _tracer->StartSpan(subcmd_str);
    auto scope = opentelemetry::trace::Tracer::WithActiveSpan(span);

    ScopedSpan ss{std::move(span), std::move(scope)};
    ss.update_attributes(config);
    ss.set_status(default_status_code);

    return ss;
  }

 private:
  [[nodiscard]] static const std::string& get_os_name() noexcept {
    static std::string s{};
    if (s.empty()) {
      try {
#ifdef HICTK_SYSTEM_NAME
        s = std::string{HICTK_SYSTEM_NAME};
        std::transform(s.begin(), s.end(), s.begin(), [](const auto c) { return std::tolower(c); });
#else
        s = "unknown";
#endif
      } catch (...) {  // NOLINT
      }
    }
    return s;
  }

  [[nodiscard]] static const std::string& get_arch() noexcept {
    static std::string s{};
    if (s.empty()) {
      try {
#ifdef HICTK_SYSTEM_PROCESSOR
        s = std::string{HICTK_SYSTEM_PROCESSOR};
        std::transform(s.begin(), s.end(), s.begin(), [](const auto c) { return std::tolower(c); });
#else
        s = "unknown";
#endif
      } catch (...) {  // NOLINT
      }
    }
    return s;
  }

  [[nodiscard]] static const std::string& get_compiler_id() noexcept {
    static std::string s{};
    if (s.empty()) {
      try {
#ifdef HICTK_CXX_COMPILER_ID
        s = std::string{HICTK_CXX_COMPILER_ID};
#else
        s = "unknown";
#endif
      } catch (...) {  // NOLINT
      }
    }
    return s;
  }

  [[nodiscard]] static const std::string& get_compiler_version() noexcept {
    static std::string s{};
    if (s.empty()) {
      try {
#ifdef HICTK_CXX_COMPILER_VERSION
        s = std::string{HICTK_CXX_COMPILER_VERSION};
#else
        s = "unknown";
#endif
      } catch (...) {  // NOLINT
      }
    }
    return s;
  }

  [[nodiscard]] static const std::string& get_build_type() noexcept {
    static std::string s{};
    if (s.empty()) {
      try {
#ifdef HICTK_BUILD_TYPE
        s = std::string{HICTK_BUILD_TYPE};
#else
        s = "unknown";
#endif
      } catch (...) {  // NOLINT
      }
    }
    return s;
  }

  [[nodiscard]] static bool init_local_telemetry_tracer() noexcept {
    try {
      namespace trace_sdk = opentelemetry::sdk::trace;
      namespace trace_exporter = opentelemetry::exporter::trace;
      std::shared_ptr<opentelemetry::trace::TracerProvider> provider =
          trace_sdk::TracerProviderFactory::Create(trace_sdk::SimpleSpanProcessorFactory::Create(
              trace_exporter::OStreamSpanExporterFactory::Create()));
      trace_api::Provider::SetTracerProvider(std::move(provider));
      return true;
    } catch (const std::exception& e) {
      SPDLOG_DEBUG(FMT_STRING("init_local_telemetry_tracer() failed: {}"), e.what());
    } catch (...) {
      SPDLOG_DEBUG("init_local_telemetry_tracer() failed: unknown error");
    }
    return false;
  }
  [[nodiscard]] static std::string get_exporter_otlp_endpoint() noexcept {
#ifdef HICTK_EXPORTER_OTLP_ENDPOINT
    constexpr std::string_view endpoint{HICTK_EXPORTER_OTLP_ENDPOINT};
    if (endpoint.find("/v1/traces") != 0) {
      return std::string{endpoint};
    }
    return fmt::format(FMT_STRING("{}/v1/traces"), endpoint);
#endif

    return "http://localhost:4318/v1/traces";
  }

  [[nodiscard]] static std::pair<std::string, std::string> get_otlp_http_header() noexcept {
#ifdef HICTK_EXPORTER_OTLP_ENDPOINT
    // NOLINTNEXTLINE(*-mt-unsafe)
    const auto* ptr = std::getenv("HICTK_EXPORTER_OTLP_HEADERS");
    if (!ptr) {
      return {};
    }

    const std::string_view header{ptr};
    const auto sep = header.find('=');
    if (sep == std::string_view::npos) {
      return {};
    }

    return std::make_pair(std::string{header.substr(0, sep)}, std::string{header.substr(sep + 1)});
#else
    return {};
#endif
  }

  [[nodiscard]] static auto get_otlp_endpoint_data() noexcept {
    struct Result {
      std::string url{};
      std::string header_name{};
      std::string header_data{};
    };

    auto [name, data] = get_otlp_http_header();
    if (name.empty()) {
      return Result{"http://localhost:4318/v1/traces", "", ""};
    }

    return Result{get_exporter_otlp_endpoint(), std::move(name), std::move(data)};
  }

  [[nodiscard]] static bool init_remote_telemetry_tracer() noexcept {
    try {
      namespace trace_sdk = opentelemetry::sdk::trace;
      namespace otlp = opentelemetry::exporter::otlp;
      namespace resource = opentelemetry::sdk::resource;

      static const std::string version{hictk::config::version::str()};

      const resource::ResourceAttributes resource_attributes{
          {"service.name", "hictk"},
          {"service.version", version},
          {"build.type", get_build_type()},
          {"build.compiler-id", get_compiler_id()},
          {"build.compiler-version", get_compiler_version()},
          {"os.type", get_os_name()},
          {"os.arch", get_arch()},
      };

      otlp::OtlpHttpExporterOptions opts{};
      const auto [url, header_name, header_data] = get_otlp_endpoint_data();

      opts.url = url;
      opts.compression = "gzip";
      if (!header_name.empty()) {
        opts.http_headers.emplace(header_name, header_data);
      }

      std::shared_ptr<trace_api::TracerProvider> provider =
          trace_sdk::TracerProviderFactory::Create(
              trace_sdk::BatchSpanProcessorFactory::Create(
                  otlp::OtlpHttpExporterFactory::Create(opts), {}),
              resource::Resource::Create(resource_attributes));
      trace_api::Provider::SetTracerProvider(std::move(provider));

      return true;
    } catch (const std::exception& e) {
      SPDLOG_DEBUG(FMT_STRING("init_remote_telemetry_tracer() failed: {}"), e.what());
    } catch (...) {
      SPDLOG_DEBUG("init_remote_telemetry_tracer() failed: unknown error");
    }
    return false;
  }

  static void init_opentelemetry_logger_once() noexcept {
    try {
      static std::once_flag flag;
      std::call_once(flag, []() {
        namespace otel_log = opentelemetry::sdk::common::internal_log;
        otel_log::GlobalLogHandler::SetLogHandler(
            std::dynamic_pointer_cast<opentelemetry::sdk::common::internal_log::LogHandler>(
                std::make_shared<OpenTelemetryLogHandler>()));

        using LogLevel = opentelemetry::sdk::common::internal_log::LogLevel;
        if constexpr (SPDLOG_ACTIVE_LEVEL < SPDLOG_LEVEL_INFO) {
          otel_log::GlobalLogHandler::SetLogLevel(LogLevel::Debug);
        } else {
          otel_log::GlobalLogHandler::SetLogLevel(LogLevel::None);
        }
      });
    } catch (...) {  // NOLINT
    }
  }

  static void register_reset_opentelemetry_tracer_at_exit() noexcept {
    try {
      auto fx = +[]() noexcept {
        trace_api::Provider::SetTracerProvider(
            std::shared_ptr<opentelemetry::trace::TracerProvider>{});
      };

      static std::once_flag flag1;
      std::call_once(flag1, &std::atexit, fx);
#ifndef _LIBCPP_VERSION
      static std::once_flag flag2;
      std::call_once(flag2, &std::at_quick_exit, fx);
#endif
    } catch (...) {  // NOLINT
    }
  }
};

#else

namespace hictk::tools {

struct ScopedSpan {
  template <typename T>
  void update_attributes([[maybe_unused]] const T& c) const noexcept {}
  template <typename T>
  void set_status([[maybe_unused]] T s) const noexcept {}
};

struct Tracer {
  enum class StatusCode : std::uint_fast8_t { kUnset, kOk, kError };

  Tracer() noexcept {
    SPDLOG_DEBUG(
        "hictk was compiled with HICTK_ENABLE_TELEMETRY=OFF: no telemetry information will be "
        "collected.");
  }
  template <typename Config>
  [[nodiscard]] std::optional<ScopedSpan> get_scoped_span(
      [[maybe_unused]] Cli::subcommand subcmd, [[maybe_unused]] const Config& config,
      [[maybe_unused]] StatusCode s = StatusCode::kError) const noexcept {
    return {};
  }
};

#endif
}
