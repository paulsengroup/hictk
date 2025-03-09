# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED)

add_library(hictk_internal_file INTERFACE)
add_library(hictk::file ALIAS hictk_internal_file)

target_sources(
  hictk_internal_file
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_file
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
target_link_libraries(
  hictk_internal_file
  INTERFACE
    hictk::balancing
    hictk::cooler
    hictk::hic
    hictk::reference
)

target_link_system_libraries(hictk_internal_file INTERFACE fmt::fmt-header-only)
