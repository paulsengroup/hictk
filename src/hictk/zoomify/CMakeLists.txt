# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)
find_package(HighFive REQUIRED QUIET)
find_package(spdlog REQUIRED QUIET)
find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_zoomify STATIC)

target_sources(
  hictk_internal_tools_zoomify
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/cool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/hic.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/mcool.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/zoomify.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/zoomify.hpp"
)

target_include_directories(hictk_internal_tools_zoomify PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_zoomify
  PRIVATE
    hictk_project_options
    hictk_project_warnings
    hictk_internal_tools_config
    hictk_internal_tools_definitions
    hictk::common
    hictk::cooler
    hictk::hic
    hictk::pixel
    hictk::tmpdir
)

target_link_system_libraries(
  hictk_internal_tools_zoomify
  PRIVATE
  fmt::fmt-header-only
  spdlog::spdlog_header_only
  std::filesystem
  HighFive
)
