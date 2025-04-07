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
    FILES
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/bin.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/bin_table.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/bin_table_fixed.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/bin_table_variable.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/bin_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/bin_table_fixed_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/bin_table_impl.hpp"
      "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/impl/bin_table_variable_impl.hpp"
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
