# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(spdlog REQUIRED QUIET)

set(HICTK_TELEMETRY_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated/include/")

if(NOT WIN32)
  file(LOCK "${HICTK_TELEMETRY_INCLUDE_DIR}" DIRECTORY GUARD FILE)
endif()

if(HICTK_ENABLE_TELEMETRY)
  message(STATUS "Building CLI tools WITH telemetry support")
  find_package(opentelemetry-cpp REQUIRED QUIET)
  find_package(fmt REQUIRED QUIET)
  find_package(phmap REQUIRED QUIET)
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/telemetry.hpp.in"
    "${HICTK_TELEMETRY_INCLUDE_DIR}/hictk/tools/telemetry.hpp"
    COPYONLY
  )
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/telemetry_impl.hpp.in"
    "${HICTK_TELEMETRY_INCLUDE_DIR}/hictk/tools/impl/telemetry_impl.hpp"
    COPYONLY
  )
else()
  message(STATUS "Building CLI tools WITHOUT telemetry support")
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/dummy_telemetry.hpp.in"
    "${HICTK_TELEMETRY_INCLUDE_DIR}/hictk/tools/telemetry.hpp"
    COPYONLY
  )
endif()

add_library(hictk_internal_tools_telemetry STATIC)

target_sources(hictk_internal_tools_telemetry PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/telemetry.cpp")

target_include_directories(hictk_internal_tools_telemetry PUBLIC "${HICTK_TELEMETRY_INCLUDE_DIR}")

target_compile_definitions(
  hictk_internal_tools_telemetry
  PRIVATE
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_ENABLE_TELEMETRY=1>"
    "$<$<BOOL:${HICTK_EXPORTER_OTLP_ENDPOINT}>:HICTK_EXPORTER_OTLP_ENDPOINT=\"${HICTK_EXPORTER_OTLP_ENDPOINT}\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_SYSTEM_NAME=\"${CMAKE_SYSTEM_NAME}\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_SYSTEM_PROCESSOR=\"${CMAKE_SYSTEM_PROCESSOR}\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_BUILD_TYPE=\"$<CONFIG>\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_CXX_COMPILER_ID=\"$<CXX_COMPILER_ID>\">"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:HICTK_CXX_COMPILER_VERSION=\"$<CXX_COMPILER_VERSION>\">"
)

target_link_libraries(
  hictk_internal_tools_telemetry
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_definitions
  PUBLIC
    hictk_internal_tools_config
    hictk_internal_tools_cli
)

target_link_system_libraries(
  hictk_internal_tools_telemetry
  PRIVATE
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:fmt::fmt-header-only>"
    "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:phmap>"
  PUBLIC
    spdlog::spdlog_header_only
)

if(TARGET opentelemetry-cpp::opentelemetry-cpp)
  target_link_libraries(
    hictk_internal_tools_telemetry
    PUBLIC
      "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:opentelemetry-cpp::opentelemetry-cpp>"
  )
else()
  target_link_system_libraries(
    hictk_internal_tools_telemetry
    PUBLIC
      "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:opentelemetry-cpp::otlp_http_exporter>"
      "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:opentelemetry-cpp::ostream_span_exporter>"
      "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:opentelemetry-cpp::sdk>"
      "$<$<BOOL:${HICTK_ENABLE_TELEMETRY}>:opentelemetry-cpp::trace>"
  )
endif()
