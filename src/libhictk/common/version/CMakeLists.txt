# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)

add_library(hictk_internal_common_version INTERFACE)
add_library(hictk::version ALIAS hictk_internal_common_version)
add_dependencies(hictk_internal_common_version _hictk_check_git)

target_sources(
  hictk_internal_common_version
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_common_version
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(hictk_internal_common_version INTERFACE fmt::fmt-header-only)
