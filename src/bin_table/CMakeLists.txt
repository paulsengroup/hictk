# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_bin_table INTERFACE)
add_library(hictk::bin_table ALIAS hictk_bin_table)
target_sources(hictk_bin_table INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/bin_table_impl.hpp)

target_include_directories(hictk_bin_table INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

target_link_libraries(
  hictk_bin_table
  INTERFACE hictk_project_warnings
            hictk_project_options
            hictk::chromosome
            hictk::common
            hictk::genomic_interval
            hictk::reference)

install(TARGETS hictk_bin_table INCLUDES)
