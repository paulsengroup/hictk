# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)

add_library(hictk_internal_formatting INTERFACE)
add_library(hictk::format ALIAS hictk_internal_formatting)

target_sources(
  hictk_internal_formatting
  INTERFACE
    FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/fmt.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/fmt/bin.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/fmt/chromosome.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/fmt/common.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/fmt/genomic_interval.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/fmt/pixel.hpp"
)

target_include_directories(
  hictk_internal_formatting
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_formatting
  INTERFACE
    hictk::bin_table
    hictk::chromosome
    hictk::genomic_interval
    hictk::pixel
    hictk::string
)

target_link_system_libraries(hictk_internal_formatting INTERFACE fmt::fmt-header-only)
