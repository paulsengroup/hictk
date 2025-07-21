# Copyright (C) 2025 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

add_library(hictk_internal_common_type_pretty_printer INTERFACE)
add_library(hictk::type_pretty_printer ALIAS hictk_internal_common_type_pretty_printer)

target_sources(
  hictk_internal_common_type_pretty_printer
  INTERFACE
    FILE_SET HEADERS
      BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
      FILES "${CMAKE_CURRENT_SOURCE_DIR}/include/hictk/type_pretty_printer.hpp"
)

target_include_directories(
  hictk_internal_common_type_pretty_printer
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:include>"
)
