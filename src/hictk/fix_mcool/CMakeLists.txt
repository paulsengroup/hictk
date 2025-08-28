# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(HighFive REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_fix_mcool STATIC)

target_sources(hictk_internal_tools_fix_mcool PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/fix_mcool.cpp")

target_include_directories(hictk_internal_tools_fix_mcool PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_fix_mcool
  PRIVATE
    hictk::project_options
    hictk::project_warnings
    hictk_internal_tools_balance
    hictk_internal_tools_validate
    hictk_internal_tools_zoomify
    hictk::cooler
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_fix_mcool
  PRIVATE
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    std::filesystem
    HighFive
)
