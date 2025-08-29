# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(HighFive REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)
find_package(Threads REQUIRED QUIET)

add_library(hictk_internal_tools_balance STATIC)

target_sources(
  hictk_internal_tools_balance
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/balance_ice.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/balance_scale.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/balance_vc.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/common.hpp"
)

target_include_directories(hictk_internal_tools_balance PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_balance
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_definitions
    hictk::balancing
    hictk::cooler
    hictk::hic
    hictk::file
    hictk::tmpdir
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_balance
  PRIVATE
    "HighFive$<$<VERSION_GREATER_EQUAL:${HighFive_VERSION},3>:::HighFive>"
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    std::filesystem
  PUBLIC
    Threads::Threads
)
