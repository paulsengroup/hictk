# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

find_package(fmt REQUIRED QUIET)

add_library(hictk_internal_common_genomic_units INTERFACE)
add_library(hictk::genomic_units ALIAS hictk_internal_common_genomic_units)

target_sources(
  hictk_internal_common_genomic_units
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_common_genomic_units
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_common_genomic_units
  INTERFACE
    hictk::common
    hictk::numeric
    hictk::type_pretty_printer
)

target_link_system_libraries(
  hictk_internal_common_genomic_units
  INTERFACE
        fmt::fmt-header-only
)
