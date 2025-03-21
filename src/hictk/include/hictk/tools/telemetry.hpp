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
#include <cstdlib>
#include <exception>
#include <memory>
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
  static void update_tracer_attributes([[maybe_unused]] SpanPtr& span,
                                       [[maybe_unused]] const std::monostate& c) {}
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
  inline static std::unique_ptr<Tracer> _instance{};

  Tracer() noexcept {
    if (!should_collect_telemetry()) {
      SPDLOG_DEBUG(
          "HICTK_NO_TELEMETRY found in environment variable list: no telemetry information will be "
          "collected.");
      return;
    }

    if (init_remote_telemetry_tracer()) {
      try {
        init_opentelemetry_logger();
        _provider = opentelemetry::trace::Provider::GetTracerProvider();
        _tracer = _provider->GetTracer("hictk");
      } catch (...) {  // NOLINT
      }
    }
  }

 public:
  using StatusCode = opentelemetry::trace::StatusCode;

  Tracer(const Tracer&) = delete;
  Tracer(Tracer&&) noexcept = delete;

  ~Tracer() noexcept {
    if (!!_tracer) {
      // NOLINTNEXTLINE(*-avoid-magic-numbers)
      _tracer->ForceFlush(std::chrono::milliseconds{500});
    }
    trace_api::Provider::SetTracerProvider(std::shared_ptr<opentelemetry::trace::TracerProvider>{});
  }

  Tracer& operator=(const Tracer&) = delete;
  Tracer& operator=(Tracer&&) noexcept = delete;

  [[nodiscard]] static Tracer* instance() noexcept {
    if (!_instance) {
      try {
        _instance = std::unique_ptr<Tracer>(new Tracer);
      } catch (...) {  // NOLINT
        return nullptr;
      }
    }
    return _instance.get();
  }

  static void tear_down_instance() noexcept { _instance.reset(); }

  template <typename Config>
  [[nodiscard]] auto get_scoped_span(Cli::subcommand subcmd, const Config& config,
                                     StatusCode default_status_code = StatusCode::kError) noexcept
      -> std::optional<ScopedSpan<SpanT, ScopeT>> {
    if (!_tracer) {
      return {};
    }

    std::string subcmd_str{Cli::subcommand_to_str(subcmd)};  // NOLINT(*-const-correctness)
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

  [[nodiscard]] static bool should_collect_telemetry() noexcept {
    // NOLINTNEXTLINE(*-implicit-bool-conversion, *-mt-unsafe)
    return !std::getenv("HICTK_NO_TELEMETRY");
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

      auto x1 = trace_exporter::OStreamSpanExporterFactory::Create();
      auto x2 = trace_sdk::SimpleSpanProcessorFactory::Create(std::move(x1));
      std::shared_ptr<opentelemetry::trace::TracerProvider> provider =
          trace_sdk::TracerProviderFactory::Create(std::move(x2));
      trace_api::Provider::SetTracerProvider(std::move(provider));
      return true;
    } catch (const std::exception& e) {
      SPDLOG_DEBUG(FMT_STRING("init_local_telemetry_tracer() failed: {}"), e.what());
    } catch (...) {
      SPDLOG_DEBUG("init_local_telemetry_tracer() failed: unknown error");
    }
    return false;
  }

  // NOLINTNEXTLINE(bugprone-exception-escape)
  [[nodiscard]] static std::string get_exporter_otlp_endpoint() noexcept {
#ifdef HICTK_EXPORTER_OTLP_ENDPOINT
    constexpr std::string_view endpoint{HICTK_EXPORTER_OTLP_ENDPOINT};
    if (endpoint.find("/v1/traces") != 0) {
      return std::string{endpoint};
    }
    return fmt::format(FMT_STRING("{}/v1/traces"), endpoint);
#else
    return "";
#endif
  }

  [[nodiscard]] static auto generate_http_exporter_opts() noexcept {
    namespace otlp = opentelemetry::exporter::otlp;
    otlp::OtlpHttpExporterOptions opts{};
    opts.url = get_exporter_otlp_endpoint();
    if (opts.url.empty()) {
      return std::unique_ptr<opentelemetry::sdk::trace::SpanExporter>{};
    }
    opts.compression = "gzip";
    opts.timeout = std::chrono::seconds{5};  // NOLINT(*-avoid-magic-numbers)
    opts.ssl_insecure_skip_verify = true;    // false;
    opts.ssl_min_tls = "1.3";
    // opts.ssl_ca_cert_string = "TODO";

    return otlp::OtlpHttpExporterFactory::Create(opts);
  }

  [[nodiscard]] static bool init_remote_telemetry_tracer() noexcept {
    try {
      namespace trace_sdk = opentelemetry::sdk::trace;
      namespace resource = opentelemetry::sdk::resource;

      auto opts = generate_http_exporter_opts();
      if (!opts) {
        return false;
      }

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

      auto pf = trace_sdk::BatchSpanProcessorFactory::Create(std::move(opts), {});
      const auto res = resource::Resource::Create(resource_attributes);
      std::shared_ptr<trace_api::TracerProvider> provider =
          trace_sdk::TracerProviderFactory::Create(std::move(pf), res);
      trace_api::Provider::SetTracerProvider(std::move(provider));

      return true;
    } catch (const std::exception& e) {
      SPDLOG_DEBUG(FMT_STRING("init_remote_telemetry_tracer() failed: {}"), e.what());
    } catch (...) {
      SPDLOG_DEBUG("init_remote_telemetry_tracer() failed: unknown error");
    }
    return false;
  }

  static void init_opentelemetry_logger() noexcept {
    try {
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
  [[nodiscard]] static Tracer* instance() noexcept { return nullptr; }
  template <typename Config>
  [[nodiscard]] std::optional<ScopedSpan> get_scoped_span(
      [[maybe_unused]] Cli::subcommand subcmd, [[maybe_unused]] const Config& config,
      [[maybe_unused]] StatusCode s = StatusCode::kError) const noexcept {
    return {};
  }
  static void tear_down_instance() noexcept {}
  [[nodiscard]] static bool should_collect_telemetry() noexcept { return false; }
};

#endif
}
