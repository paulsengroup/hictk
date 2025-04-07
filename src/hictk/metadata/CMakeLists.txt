# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(nlohmann_json REQUIRED QUIET)
find_package(phmap REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_metadata STATIC)

target_sources(
  hictk_internal_tools_metadata
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/common.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/cool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mcool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/metadata.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/scool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/metadata.hpp"
)

target_include_directories(hictk_internal_tools_metadata PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_metadata
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_definitions
    hictk_internal_tools_file_attributes_formatting
    hictk::bin_table
    hictk::common
    hictk::cooler
    hictk::hic
  PUBLIC
    hictk_internal_tools_config
)

target_link_system_libraries(
  hictk_internal_tools_metadata
  PRIVATE
    fmt::fmt-header-only
    nlohmann_json::nlohmann_json
    std::filesystem
    phmap
)
