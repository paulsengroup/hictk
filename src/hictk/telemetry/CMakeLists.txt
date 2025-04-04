# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(phmap REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

if(HICTK_ENABLE_TELEMETRY)
  message(STATUS "Building with telemetry support")
  find_package(opentelemetry-cpp REQUIRED QUIET)
  find_package(fmt REQUIRED QUIET)
  find_package(phmap REQUIRED QUIET)
else()
  message(STATUS "Building without telemetry support")
endif()

add_library(hictk_internal_tools_telemetry INTERFACE)

target_include_directories(hictk_internal_tools_telemetry INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_compile_definitions(
  hictk_internal_tools_telemetry
  INTERFACE
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_ENABLE_TELEMETRY>"
    "$<$<BOOL:${HICTK_EXPORTER_OTLP_ENDPOINT}>:HICTK_EXPORTER_OTLP_ENDPOINT=\"${HICTK_EXPORTER_OTLP_ENDPOINT}\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_SYSTEM_NAME=\"${CMAKE_SYSTEM_NAME}\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_SYSTEM_PROCESSOR=\"${CMAKE_SYSTEM_PROCESSOR}\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_BUILD_TYPE=\"$<CONFIG>\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_CXX_COMPILER_ID=\"$<CXX_COMPILER_ID>\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_CXX_COMPILER_VERSION=\"$<CXX_COMPILER_VERSION>\">"
)

target_link_libraries(
  hictk_internal_tools_telemetry
  INTERFACE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_config
    hictk_internal_tools_definitions
    hictk_internal_tools_cli
)

target_link_system_libraries(
  hictk_internal_tools_telemetry
  INTERFACE
  "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:fmt::fmt-header-only>"
  "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:opentelemetry-cpp::opentelemetry-cpp>"
  spdlog::spdlog_header_only
  std::filesystem
  phmap
)
