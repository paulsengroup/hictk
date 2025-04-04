# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_bin_table INTERFACE)
add_library(hictk::bin_table ALIAS hictk_internal_bin_table)
target_sources(
  hictk_internal_bin_table
  INTERFACE
    FILE_SET HEADERS
    BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_include_directories(
  hictk_internal_bin_table
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)

target_link_libraries(
  hictk_internal_bin_table
  INTERFACE
    hictk::chromosome
    hictk::common
    hictk::genomic_interval
    hictk::reference
    hictk::type_traits
)
