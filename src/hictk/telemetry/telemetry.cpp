// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef HICTK_ENABLE_TELEMETRY

// clang-format off
#include "hictk/tools/telemetry.hpp"
// clang-format on

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
#include <spdlog/spdlog.h>

#include <cassert>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

#include "hictk/tools/cli.hpp"
#include "hictk/tools/config.hpp"

namespace hictk::tools {

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

Tracer::Tracer() noexcept {
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

Tracer::~Tracer() noexcept {
  if (!!_tracer) {
    // NOLINTNEXTLINE(*-avoid-magic-numbers)
    _tracer->ForceFlush(std::chrono::milliseconds{500});
  }
  trace_api::Provider::SetTracerProvider(std::shared_ptr<opentelemetry::trace::TracerProvider>{});
}

Tracer* Tracer::instance() noexcept {
  if (!_instance) {
    try {
      _instance = std::unique_ptr<Tracer>(new Tracer);
    } catch (...) {  // NOLINT
      return nullptr;
    }
  }
  return _instance.get();
}

void Tracer::tear_down_instance() noexcept { _instance.reset(); }

bool Tracer::should_collect_telemetry() noexcept {
  // NOLINTNEXTLINE(*-implicit-bool-conversion, *-mt-unsafe)
  return !std::getenv("HICTK_NO_TELEMETRY");
}

const std::string& Tracer::get_os_name() noexcept {
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

const std::string& Tracer::get_arch() noexcept {
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

const std::string& Tracer::get_compiler_id() noexcept {
  static std::string s{};
  if (s.empty()) {
    try {
#ifdef HICTK_CXX_COMPILER_ID
      s = std::string{HICTK_CXX_COMPILER_ID};
      std::transform(s.begin(), s.end(), s.begin(), [](const auto c) { return std::tolower(c); });
#else
      s = "unknown";
#endif
    } catch (...) {  // NOLINT
    }
  }
  return s;
}

const std::string& Tracer::get_compiler_version() noexcept {
  static std::string s{};
  if (s.empty()) {
    try {
#ifdef HICTK_CXX_COMPILER_VERSION
      s = std::string{HICTK_CXX_COMPILER_VERSION};
      std::transform(s.begin(), s.end(), s.begin(), [](const auto c) { return std::tolower(c); });
#else
      s = "unknown";
#endif
    } catch (...) {  // NOLINT
    }
  }
  return s;
}

const std::string& Tracer::get_build_type() noexcept {
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

bool Tracer::init_local_telemetry_tracer() noexcept {
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
std::string Tracer::get_exporter_otlp_endpoint() noexcept {
#ifdef HICTK_EXPORTER_OTLP_ENDPOINT
  constexpr std::string_view endpoint{HICTK_EXPORTER_OTLP_ENDPOINT};
  if (endpoint.find("/v1/traces") != std::string_view::npos) {
    return std::string{endpoint};
  }
  return fmt::format(FMT_STRING("{}/v1/traces"), endpoint);
#else
  return "";
#endif
}

auto Tracer::generate_http_exporter_opts() noexcept {
  namespace otlp = opentelemetry::exporter::otlp;
  otlp::OtlpHttpExporterOptions opts{};
  opts.url = get_exporter_otlp_endpoint();
  if (opts.url.empty()) {
    return std::unique_ptr<opentelemetry::sdk::trace::SpanExporter>{};
  }
  opts.compression = "gzip";
  opts.timeout = std::chrono::seconds{5};  // NOLINT(*-avoid-magic-numbers)
  opts.ssl_insecure_skip_verify = false;
  opts.ssl_min_tls = "1.2";

  return otlp::OtlpHttpExporterFactory::Create(opts);
}

bool Tracer::init_remote_telemetry_tracer() noexcept {
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

void Tracer::init_opentelemetry_logger() noexcept {
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
}  // namespace hictk::tools

#endif
