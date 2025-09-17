// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#ifdef HICTK_ENABLE_TELEMETRY

// clang-format off
#include "hictk/tools/telemetry.hpp"
// clang-format on

#include <fmt/format.h>
#include <fmt/ranges.h>
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
#include <spdlog/spdlog.h>

#include <boost/hash2/sha3.hpp>
#include <cassert>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <memory>
#include <nonstd/span.hpp>
#include <string>
#include <utility>

#include "hictk/tools/build_options.hpp"
#include "hictk/tools/cli.hpp"

namespace hictk::tools {

namespace internal {

std::string hash_argv(nonstd::span<const char*> argv) noexcept {
  try {
    if (argv.empty()) {
      return {};
    }

    boost::hash2::sha3_256 hasher;
    for (const char* arg : argv) {
      hasher.update(arg, std::strlen(arg));
    }

    const auto digest = hasher.result();
    return fmt::format(FMT_STRING("{:x}"), fmt::join(digest.begin(), digest.end(), ""));

  } catch (const std::exception& e) {
    SPDLOG_WARN(FMT_STRING("failed to hash command line arguments: {}"), e.what());
  } catch (...) {
    SPDLOG_WARN("failed to hash command line arguments: unknown error");
  }

  return {};
}
}  // namespace internal

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
    // NOLINTBEGIN
    try {
      _instance = new Tracer{};
    } catch (...) {
      return nullptr;
    }
    // NOLINTEND
  }
  return _instance;
}

void Tracer::tear_down_instance() noexcept {
  delete _instance;  // NOLINT
  _instance = nullptr;
}

bool Tracer::should_collect_telemetry() noexcept {
  // NOLINTNEXTLINE(*-implicit-bool-conversion, *-mt-unsafe)
  return !std::getenv("HICTK_NO_TELEMETRY");
}

[[nodiscard]] static opentelemetry::sdk::resource::ResourceAttributes
generate_resource_attributes() {
  auto try_get_build_attr =
      [build = get_build_options_json()](std::string_view key) noexcept -> std::string {
    try {
      return build.at(key).get<std::string>();
    } catch (...) {  // NOLINT
      return "unknown";
    }
  };

  opentelemetry::sdk::resource::ResourceAttributes attrs{
      {"service.name", "hictk"},
      {"service.version", std::string{config::version::str()}},
      {"host.arch", try_get_build_attr("arch")},
      {"build.compiler.name", try_get_build_attr("compiler_name")},
      {"build.compiler.version", try_get_build_attr("compiler_version")},
      {"build.type", try_get_build_attr("build_type")},
      {"os.type", try_get_build_attr("os_name")},
      {"os.version", try_get_build_attr("os_version")},
  };

  const auto deps = get_dependency_versions_json();
  for (const auto& [key, value] : deps.items()) {
    std::string name = key;
    std::transform(name.begin(), name.end(), name.begin(),
                   [](const auto c) { return std::tolower(c); });
    attrs.emplace(fmt::format(FMT_STRING("build.dependencies.{}.version"), name),
                  value.get<std::string>());
  }

  return attrs;
}

bool Tracer::init_local_telemetry_tracer() noexcept {
  try {
    namespace trace_sdk = opentelemetry::sdk::trace;
    namespace trace_exporter = opentelemetry::exporter::trace;
    namespace resource = opentelemetry::sdk::resource;

    auto x1 = trace_exporter::OStreamSpanExporterFactory::Create();
    auto x2 = trace_sdk::SimpleSpanProcessorFactory::Create(std::move(x1));
    const auto res = resource::Resource::Create(generate_resource_attributes());
    std::shared_ptr<opentelemetry::trace::TracerProvider> provider =
        trace_sdk::TracerProviderFactory::Create(std::move(x2), res);
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

    auto pf = trace_sdk::BatchSpanProcessorFactory::Create(std::move(opts), {});
    const auto res = resource::Resource::Create(generate_resource_attributes());
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
