// Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include "hictk/tools/dependency_metadata.hpp"

#include <H5public.h>
#include <archive.h>
#include <fmt/format.h>
#include <spdlog/spdlog.h>
#include <zstd.h>

#include <exception>
#include <nlohmann/json.hpp>
#include <string>
#include <string_view>

namespace hictk::tools {

[[nodiscard]] static std::string get_hdf5_library_version() noexcept {
  try {
    unsigned major{};
    unsigned minor{};
    unsigned patch{};

    const auto status = H5get_libversion(&major, &minor, &patch);
    if (status == 0) {
      return fmt::format(FMT_STRING("{}.{}.{}"), major, minor, patch);
    }

  } catch (...) {  // NOLINT
  }
  return "unknown";
}

[[nodiscard]] static std::string_view get_libarchive_library_version() noexcept {
  constexpr std::string_view prefix{"libarchive "};
  std::string_view version{archive_version_string()};
  version.remove_prefix(prefix.size());

  return version;
}

nlohmann::json get_dependency_versions_json() noexcept {
  try {
    nlohmann::json deps;

    // compile-time deps
    deps.emplace("boost", HICTK_BOOST_VERSION);
    deps.emplace("bshoshany-thread-pool", HICTK_BSHOSHANY_THREAD_POOL_VERSION);
    deps.emplace("CLI11", HICTK_CLI11_VERSION);
    deps.emplace("concurrentqueue", HICTK_CONCURRENTQUEUE_VERSION);
    deps.emplace("fast_float", HICTK_FASTFLOAT_VERSION);
    deps.emplace("fmt", HICTK_FMT_VERSION);
    deps.emplace("HighFive", HICTK_HIGHFIVE_VERSION);
    deps.emplace("nlohmann_json", HICTK_NLOHMANN_JSON_VERSION);
    deps.emplace("parallel-hashmap", HICTK_PHMAP_VERSION);
    deps.emplace("readerwriterqueue", HICTK_READERWRITERQUEUE_VERSION);
    deps.emplace("span-lite", HICTK_SPAN_LITE_VERSION);
    deps.emplace("tomlplusplus", HICTK_TOMLPLUSPLUS_VERSION);

    // runtime-deps
    deps.emplace("HDF5", get_hdf5_library_version());
    deps.emplace("LibArchive", get_libarchive_library_version());
    deps.emplace("libdeflate", HICTK_LIBDEFLATE_VERSION);
    deps.emplace("opentelemetry-cpp", HICTK_OPENTELEMETRY_CPP_VERSION);
    deps.emplace("zstd", ZSTD_versionString());

    return deps;

  } catch (const std::exception& e) {
    SPDLOG_WARN(FMT_STRING("failed to collect dependency versions: {}"), e.what());
  } catch (...) {
    SPDLOG_WARN("failed to collect dependency versions: unknown error");
  }

  return {};
}

std::string get_dependency_versions(bool pretty) noexcept {
  try {
    if (pretty) {
      return get_dependency_versions_json().dump(2, ' ');
    }
    return get_dependency_versions_json().dump();
  } catch (const std::exception& e) {
    SPDLOG_WARN(FMT_STRING("failed to collect dependency versions: {}"), e.what());
  } catch (...) {
    SPDLOG_WARN("failed to collect dependency versions: unknown error");
  }

  return "{}";
}

}  // namespace hictk::tools
