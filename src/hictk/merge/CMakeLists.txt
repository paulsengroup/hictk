# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_merge STATIC)

target_sources(
  hictk_internal_tools_merge
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/cool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/merge.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/merge.hpp"
)

target_include_directories(hictk_internal_tools_merge PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_merge
  PRIVATE
    hictk::project_options
    hictk::project_warnings
    hictk::file
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_merge
  PRIVATE
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    std::filesystem
)
