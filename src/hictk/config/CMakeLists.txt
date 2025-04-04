# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(Filesystem REQUIRED QUIET)

add_library(hictk_internal_tools_config INTERFACE)

target_include_directories(hictk_internal_tools_config INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")

target_link_libraries(
  hictk_internal_tools_config
  INTERFACE
    hictk_project_options
    hictk_project_warnings
    hictk::common
    hictk::balancing
    hictk::hic
)

target_link_system_libraries(
        hictk_internal_tools_config INTERFACE std::filesystem
)
