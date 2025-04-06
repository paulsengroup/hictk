# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_validate STATIC)

target_sources(
  hictk_internal_tools_validate
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/cool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mcool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/scool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/validate.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/validate.hpp"
)

target_include_directories(hictk_internal_tools_validate PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_validate
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_definitions
    hictk_internal_tools_file_attributes_formatting
    hictk::chromosome
    hictk::cooler
    hictk::hic
    hictk::string
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_validate
  PRIVATE
    fmt::fmt-header-only
    spdlog::spdlog_header_only
    std::filesystem
)
