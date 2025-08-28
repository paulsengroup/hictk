# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)

add_library(hictk_internal_tools_runner STATIC)

target_sources(hictk_internal_tools_runner PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/runner.cpp")

target_include_directories(hictk_internal_tools_runner PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_runner
  PRIVATE
    hictk::project_options
    hictk::project_warnings
    hictk_internal_tools_balance
    hictk_internal_tools_convert
    hictk_internal_tools_dump
    hictk_internal_tools_fix_mcool
    hictk_internal_tools_load
    hictk_internal_tools_merge
    hictk_internal_tools_metadata
    hictk_internal_tools_rename_chromosomes
    hictk_internal_tools_telemetry
    hictk_internal_tools_validate
    hictk_internal_tools_zoomify
  PUBLIC
    hictk_internal_tools_cli
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_runner
  PRIVATE
    fmt::fmt-header-only
    spdlog::spdlog_header_only
)
