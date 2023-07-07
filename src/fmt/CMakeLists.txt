# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_fmt INTERFACE)
add_library(hictk::fmt ALIAS hictk_fmt)

target_include_directories(hictk_fmt INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/)

target_link_libraries(
  hictk_fmt
  INTERFACE hictk_project_warnings
            hictk_project_options
            hictk::bin_table
            hictk::chromosome
            hictk::common
            hictk::genomic_interval
            hictk::pixel)

target_link_system_libraries(hictk_fmt INTERFACE _hictk_fmt_tgt)

install(TARGETS hictk_fmt INCLUDES)
