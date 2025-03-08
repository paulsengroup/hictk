# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(FMT REQUIRED)

add_library(hictk_formatting INTERFACE)
add_library(hictk::format ALIAS hictk_formatting)

target_sources(
  hictk_formatting
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_formatting
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_formatting
  INTERFACE
    hictk::bin_table
    hictk::chromosome
    hictk::common
    hictk::genomic_interval
    hictk::pixel
)

target_link_system_libraries(hictk_formatting INTERFACE fmt::fmt-header-only)
